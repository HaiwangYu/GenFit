#include <ConstField.h>
#include <Exception.h>
#include <FieldManager.h>
#include <KalmanFitterRefTrack.h>
#include <StateOnPlane.h>
#include <Track.h>
#include <TrackPoint.h>

#include <MaterialEffects.h>
#include <RKTrackRep.h>
#include <G4eTrackRep.h>
#include <TGeoMaterialInterface.h>

#include <EventDisplay.h>

#include <HelixTrackModel.h>
#include <MeasurementCreator.h>

#include <TDatabasePDG.h>
#include <TEveManager.h>
#include <TGeoManager.h>
#include <TRandom.h>
#include <TVector3.h>
#include <vector>

#include "TDatabasePDG.h"
#include <TMath.h>

#include "PlanarMeasurement.h"
#include "DetPlane.h"
#include "SharedPlanePtr.h"
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

#define LogDEBUG    std::cout<<"DEBUG: "<<__LINE__<<"\n"

int main() {

	gRandom->SetSeed(14);

	// init MeasurementCreator
	genfit::MeasurementCreator measurementCreator;

	// init geometry and mag. field
	new TGeoManager("Geometry", "Geane geometry");
	TGeoManager::Import("genfitGeom.root");
	genfit::FieldManager::getInstance()->init(
			new genfit::ConstField(40., 0., 0.)); // kGauss
	genfit::MaterialEffects::getInstance()->init(
			new genfit::TGeoMaterialInterface());

	// init event display
	genfit::EventDisplay* display = genfit::EventDisplay::getInstance();

	// init fitter
	genfit::AbsKalmanFitter* fitter = new genfit::KalmanFitterRefTrack();

	// main loop
	for (unsigned int iEvent = 0; iEvent < 100; ++iEvent) {

		// true start values
		TVector3 pos(0, 0, 0); //cm
		TVector3 mom(0, 0, 1); //GeV
//    mom.SetPhi(gRandom->Uniform(0.,2*TMath::Pi()));
//    mom.SetTheta(gRandom->Uniform(0.4*TMath::Pi(),0.6*TMath::Pi()));
//    mom.SetMag(gRandom->Uniform(0.2, 1.));

		// helix track model
		const int pdg = -13;               // particle pdg code
//		const double charge = TDatabasePDG::Instance()->GetParticle(pdg)->Charge() / (3.);
//    genfit::HelixTrackModel* helix = new genfit::HelixTrackModel(pos, mom, charge);
//    measurementCreator.setTrackModel(helix);

		unsigned int nMeasurements = gRandom->Uniform(5, 15);
		nMeasurements = 8;

		// smeared start values
		const bool smearPosMom = false; // init the Reps with smeared pos and mom
		const double posSmear = 0.1;     // cm
		const double momSmear = 3. / 180. * TMath::Pi();     // rad
		const double momMagSmear = 0.1;   // relative

		TVector3 posM(pos);
		TVector3 momM(mom);
		if (smearPosMom) {
			posM.SetX(gRandom->Gaus(posM.X(), posSmear));
			posM.SetY(gRandom->Gaus(posM.Y(), posSmear));
			posM.SetZ(gRandom->Gaus(posM.Z(), posSmear));

			momM.SetPhi(gRandom->Gaus(mom.Phi(), momSmear));
			momM.SetTheta(gRandom->Gaus(mom.Theta(), momSmear));
			momM.SetMag(gRandom->Gaus(mom.Mag(), momMagSmear * mom.Mag()));
		}
		// approximate covariance
		TMatrixDSym covM(6);
		//double resolution = 0.01;//original
		double resolution = 0.1;   //original
		measurementCreator.setResolution(resolution);
		for (int i = 0; i < 6; ++i)
			covM(i, i) = resolution * resolution;
//		for (int i = 3; i < 6; ++i)
//			covM(i, i) = pow(resolution / nMeasurements / sqrt(3), 2);

		// trackrep
		genfit::AbsTrackRep* rep = new genfit::RKTrackRep(pdg);

		// smeared start state
		genfit::MeasuredStateOnPlane stateSmeared(rep);
		stateSmeared.setPosMomCov(posM, momM, covM);

		// create track
		TVectorD seedState(6);
		TMatrixDSym seedCov(6);
		stateSmeared.get6DStateCov(seedState, seedCov);
		genfit::Track fitTrack(rep, seedState, seedCov);

		// create random measurement types
//		std::vector<genfit::eMeasurementType> measurementTypes;
//		for (unsigned int i = 0; i < nMeasurements; ++i)
//			measurementTypes.push_back(genfit::eMeasurementType(3));
		//measurementTypes.push_back(genfit::eMeasurementType(gRandom->Uniform(8)));

		genfit::MeasuredStateOnPlane currentState = stateSmeared;

		genfit::SharedPlanePtr detPlane = boost::make_shared<genfit::DetPlane>(
				TVector3(0, 0, 25), TVector3(0, 0, 1));
		//detPlane->Print();
		//rep->setDebugLvl(10);

		//currentState.Print();
		//rep->extrapolateToPlane(currentState, detPlane);
		//currentState.Print();

		int measurementCounter_ = 0;
		// create smeared measurements and add to track
		try {
			for (unsigned int i = 1; i < 5; ++i) {

				genfit::SharedPlanePtr plane(
						new genfit::DetPlane(TVector3(0, 0, i * 20),
								TVector3(0, 0, 1)));


				genfit::MeasuredStateOnPlane tempState = stateSmeared;
				rep->extrapolateToPlane(tempState, plane);
				TVectorD tempPosMom(6);
				TMatrixDSym tempPosMomCov(6);
				tempState.get6DStateCov(tempPosMom, tempPosMomCov);

				std::cout<<"DEBUG: "<<" pos: "
						<<tempPosMom(0)<<","<<tempPosMom(1)<<","<<tempPosMom(2)<<","
						<<tempPosMom(3)<<","<<tempPosMom(4)<<","<<tempPosMom(5)<<","
						<<"\n";

				int nDim = 2;
				TVectorD hitCoords(nDim);
				TMatrixDSym hitCov(2);

				for(int i0 = 0;i0<2;i0++)
				{
					hitCoords(i0) = tempPosMom(i0) + gRandom->Gaus(0, resolution);
					hitCov(i0,i0) = resolution * resolution;
				}
				genfit::AbsMeasurement* measurement =
						new genfit::PlanarMeasurement(hitCoords, hitCov,
								int(genfit::StripUV), measurementCounter_,
								nullptr);
				std::vector<genfit::AbsMeasurement*> measurements;
				static_cast<genfit::PlanarMeasurement*>(measurement)->setPlane(plane, measurementCounter_);
				measurements.push_back(measurement);
//				std::vector<genfit::AbsMeasurement*> measurements = measurementCreator.create(genfit::StripUV, i*20.);
				fitTrack.insertPoint(
						new genfit::TrackPoint(measurements, &fitTrack));
			}
		} catch (genfit::Exception& e) {
			std::cerr << "Exception, next track" << std::endl;
			std::cerr << e.what();
			continue;
		}

		//check
		assert(fitTrack.checkConsistency());

		// do the fit
		fitter->processTrack(&fitTrack);

		//check
		assert(fitTrack.checkConsistency());

		if (iEvent < 1000) {
			// add track to event display
			display->addEvent(&fitTrack);
		}

		fitTrack.getCardinalRep()->extrapolateToPoint(currentState,TVector3(0,0,0));
		std::cout<<"DEBUG: extrapolateToPoint(0,0,0)\n";
		currentState.Print();

	}    // end loop over events

	delete fitter;

	// open event display
	display->open();

}

