#include <ConstField.h>
#include <Exception.h>
#include <FieldManager.h>
#include <KalmanFitterRefTrack.h>
#include <StateOnPlane.h>
#include <Track.h>
#include <TrackPoint.h>

#include <MaterialEffects.h>
#include <RKTrackRep.h>
//#include <G4eTrackRep.h>
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
#include <KalmanFittedStateOnPlane.h>
//#include <AbsKalmanFitter.h>
//#include <KalmanFitter.h>
//#include <KalmanFitterRefTrack.h>
#include <KalmanFitterInfo.h>
//#include <KalmanFitStatus.h>
#include "TH1D.h"
#include <TFile.h>
#include <TCanvas.h>

#define LogDEBUG    std::cout<<"DEBUG: "<<__LINE__<<"\n"

int main(int argc, char **argv) {

	double magnetic_field = atof(argv[1]); //kGauss
	double momentum_z = atof(argv[2]); // GeV
	double init_z = atof(argv[3]); //cm
	double dest_z = atof(argv[4]); //cm
	std::cout << magnetic_field << "\n";
	double resolution_detector = 0.1;
	unsigned int nMeasurements = 4; //4 stations
	gRandom->SetSeed(14);

	// init MeasurementCreator
	//genfit::MeasurementCreator measurementCreator;

	// init geometry and mag. field
	new TGeoManager("Geometry", "Geane geometry");
	TGeoManager::Import("genfitGeom.root");
	genfit::FieldManager::getInstance()->init(
			new genfit::ConstField(magnetic_field, 0., 0.)); // kGauss
	genfit::MaterialEffects::getInstance()->init(
			new genfit::TGeoMaterialInterface());

	// init event display
	genfit::EventDisplay* display = genfit::EventDisplay::getInstance();

	// init fitter
	genfit::AbsKalmanFitter* fitter = new genfit::KalmanFitterRefTrack();

	TH1D *hmomRes = new TH1D("hmomRes", "mom res", 500,
			-20 * resolution_detector * momentum_z / nMeasurements,
			20 * resolution_detector * momentum_z / nMeasurements);
	TH1D *hupRes = new TH1D("hupRes", "u' res", 500,
			-15 * resolution_detector / nMeasurements,
			15 * resolution_detector / nMeasurements);
	TH1D *hvpRes = new TH1D("hvpRes", "v' res", 500,
			-15 * resolution_detector / nMeasurements,
			15 * resolution_detector / nMeasurements);
	TH1D *huRes = new TH1D("huRes", "u res", 500, -15 * resolution_detector,
			15 * resolution_detector);
	TH1D *hvRes = new TH1D("hvRes", "v res", 500, -15 * resolution_detector,
			15 * resolution_detector);

	TH1D *hqopPu = new TH1D("hqopPu", "q/p pull", 200, -6., 6.);
//	TH1D *pVal = new TH1D("pVal", "p-value", 100, 0., 1.00000001);
//	pVal->SetMinimum(0);
	TH1D *hupPu = new TH1D("hupPu", "u' pull", 200, -6., 6.);
	TH1D *hvpPu = new TH1D("hvpPu", "v' pull", 200, -6., 6.);
	TH1D *huPu = new TH1D("huPu", "u pull", 200, -6., 6.);
	TH1D *hvPu = new TH1D("hvPu", "v pull", 200, -6., 6.);

	TH1D *hchi2Forward = new TH1D("hchi2Forward", "Fit #chi^{2}", 100,
			0, 3.);

	// main loop
	for (unsigned int iEvent = 0; iEvent < 1000; ++iEvent) {

		// true start values
		TVector3 init_pos(0, 0, init_z); //cm
		TVector3 init_mom(0, 0, momentum_z); //GeV
//    init_mom.SetPhi(gRandom->Uniform(0.,2*TMath::Pi()));
//    init_mom.SetTheta(gRandom->Uniform(0.4*TMath::Pi(),0.6*TMath::Pi()));
//    init_mom.SetMag(gRandom->Uniform(0.2, 1.));

		// helix track model
		const int pdg = 13; //-13: mu+, 13: mu-
		const double charge =
				TDatabasePDG::Instance()->GetParticle(pdg)->Charge() / (3.);
//    genfit::HelixTrackModel* helix = new genfit::HelixTrackModel(init_pos, init_mom, charge);
//    measurementCreator.setTrackModel(helix);

		// Seed: use smeared values
		const bool smearPosMom = true; // init the Reps with smeared init_pos and init_mom
		const double posSmear = 10 * resolution_detector;     // cm
		const double momSmear = 3. / 180. * TMath::Pi();     // rad
		const double momMagSmear = 0.1;   // relative

		TVector3 seed_pos(init_pos);
		TVector3 seed_mom(init_mom);
		if (smearPosMom) {
			seed_pos.SetX(gRandom->Gaus(seed_pos.X(), posSmear));
			seed_pos.SetY(gRandom->Gaus(seed_pos.Y(), posSmear));
			seed_pos.SetZ(gRandom->Gaus(seed_pos.Z(), posSmear));

			seed_mom.SetPhi(gRandom->Gaus(init_mom.Phi(), momSmear));
			seed_mom.SetTheta(gRandom->Gaus(init_mom.Theta(), momSmear));
			seed_mom.SetMag(
					gRandom->Gaus(init_mom.Mag(),
							momMagSmear * init_mom.Mag()));
		}

		// approximate covariance
		TMatrixDSym seed_cov(6);

//		measurementCreator.setResolution(resolution_detector);

		for (int i = 0; i < 3; ++i)
			seed_cov(i, i) = resolution_detector * resolution_detector;
		for (int i = 3; i < 6; ++i)
			seed_cov(i, i) = pow(resolution_detector / nMeasurements / sqrt(3),
					2);

//		seed_cov(0, 0) = 0;
//		seed_cov(1, 1) = 0;
//		seed_cov(2, 2) = 0;
//		seed_cov(3, 3) = 0;
//		seed_cov(4, 4) = 0;
//		seed_cov(5, 5) = 0;

		// trackrep
		genfit::AbsTrackRep* rep = new genfit::RKTrackRep(pdg);
		genfit::AbsTrackRep* measurementCreatorRep = new genfit::RKTrackRep(pdg);


		// seed start state
		genfit::MeasuredStateOnPlane seedMSoP(rep);
		seedMSoP.setPosMomCov(seed_pos, seed_mom, seed_cov);
		const genfit::StateOnPlane seedSoP(seedMSoP);

		genfit::MeasuredStateOnPlane initMSoP(measurementCreatorRep);
		initMSoP.setPosMomCov(init_pos, init_mom, seed_cov);
		const genfit::StateOnPlane initSoP(initMSoP);


		// create track
		TVectorD seedState(6);
		TMatrixDSym seedCov(6);
		seedMSoP.get6DStateCov(seedState, seedCov);
		genfit::Track fitTrack(rep, seedState, seedCov);

		TVectorD initState(6);
		TMatrixDSym initCov(6);
		initMSoP.get6DStateCov(seedState, seedCov);
		genfit::Track measurementCreatorTrack(measurementCreatorRep, initState, initCov);

		genfit::MeasuredStateOnPlane currentState = seedMSoP;

		genfit::SharedPlanePtr destPlane = boost::make_shared<genfit::DetPlane>(
				TVector3(0, 0, dest_z), TVector3(0, 0, 1));
		//rep->setDebugLvl(10);

//		std::cout << "============= State at z = " << init_z
//				<< ": =============\n";
//		currentState.Print();
//		rep->extrapolateToPlane(currentState, destPlane);
//		std::cout << "============= State at z = " << dest_z
//				<< ": =============\n";
//		currentState.Print();

		int measurementCounter_ = 0;
		// create smeared measurements and add to track
		try {
			for (unsigned int i = 1; i < 5; ++i) {

				genfit::SharedPlanePtr plane(
						new genfit::DetPlane(TVector3(0, 0, i * 20),
								TVector3(0, 0, 1)));

//				genfit::MeasuredStateOnPlane tempState = seedMSoP;
//				rep->extrapolateToPlane(tempState, plane);
				genfit::MeasuredStateOnPlane tempState = initMSoP;
				measurementCreatorRep->extrapolateToPlane(tempState, plane);
				TVectorD tempPosMom(6);
				TMatrixDSym tempPosMomCov(6);
				tempState.get6DStateCov(tempPosMom, tempPosMomCov);

				std::cout << "DEBUG: " << " position (cm):" << tempPosMom(0)
						<< "," << tempPosMom(1) << "," << tempPosMom(2)
						<< "; \n" << " momentum (GeV): " << tempPosMom(3) << ","
						<< tempPosMom(4) << "," << tempPosMom(5) << "." << "\n";

				int nDim = 2;
				TVectorD hitCoords(nDim);
				TMatrixDSym hitCov(nDim);

				for (int i0 = 0; i0 < nDim; i0++) {
					hitCoords(i0) = tempPosMom(i0)
							+ gRandom->Gaus(0, resolution_detector);
					hitCov(i0, i0) = resolution_detector * resolution_detector;
				}
				genfit::AbsMeasurement* measurement =
						new genfit::PlanarMeasurement(hitCoords, hitCov, -1,
								measurementCounter_,
								nullptr);
//				measurement->Print();
				std::vector<genfit::AbsMeasurement*> measurements;
				static_cast<genfit::PlanarMeasurement*>(measurement)->setPlane(
						plane, measurementCounter_);
				measurements.push_back(measurement);

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

		if (iEvent < 10) {
			// add track to event display
			display->addEvent(&fitTrack);
		}

		// check if fit was successful
		if (!fitTrack.getFitStatus(rep)->isFitConverged()) {
			std::cout
					<< "Track could not be fitted successfully! Fit is not converged! \n";
			continue;
		}

//		fitTrack.Print();
		double chi2 = fitTrack.getFitStatus(rep)->getChi2();
		double ndf = fitTrack.getFitStatus(rep)->getNdf();
		hchi2Forward->Fill(chi2 / ndf);

//		fitTrack.getCardinalRep()->extrapolateToPoint(currentState,
//				TVector3(0, 0, 0));
//		std::cout << "DEBUG: extrapolateToPoint(0,0,0)\n";
//		fitTrack.getCardinalRep()
		rep->extrapolateToPlane(currentState, seedSoP.getPlane());
		std::cout << "DEBUG: yuhw extrapolateToPlane: \n";
//		currentState.Print();

		genfit::TrackPoint* tp = fitTrack.getPointWithMeasurementAndFitterInfo(
				0, rep);
		if (tp == NULL) {
			std::cout << "Track has no TrackPoint with fitterInfo! \n";
			continue;
		}
		genfit::KalmanFittedStateOnPlane kfsop(
				*(static_cast<genfit::KalmanFitterInfo*>(tp->getFitterInfo(rep))->getBackwardUpdate()));
		// extrapolate back to reference plane.
		try {
			rep->extrapolateToPlane(kfsop, initSoP.getPlane());
			std::cout << "DEBUG: Official extrapolateToPlane: \n";
//			kfsop.Print();
		} catch (genfit::Exception& e) {
			std::cerr << "Exception, next track" << std::endl;
			std::cerr << e.what();
			continue;
		}

		const TVectorD& referenceState = initSoP.getState();
		const TVectorD& state = kfsop.getState();
		const TMatrixDSym& cov = kfsop.getCov();

		hmomRes->Fill((-charge / state[0] - seed_mom.Mag()));
		hupRes->Fill((state[1] - referenceState[1]));
		hvpRes->Fill((state[2] - referenceState[2]));
		huRes->Fill((state[3] - referenceState[3]));
		hvRes->Fill((state[4] - referenceState[4]));

		hqopPu->Fill((-state[0] - referenceState[0]) / sqrt(cov[0][0]));
		hupPu->Fill((state[1] - referenceState[1]) / sqrt(cov[1][1]));
		hvpPu->Fill((state[2] - referenceState[2]) / sqrt(cov[2][2]));
		huPu->Fill((state[3] - referenceState[3]) / sqrt(cov[3][3]));
		hvPu->Fill((state[4] - referenceState[4]) / sqrt(cov[4][4]));

	}    // end loop over events

	TCanvas* c1 = new TCanvas();
	c1->Divide(2, 3);

	c1->cd(1);
	hmomRes->Fit("gaus");
	hmomRes->Draw();

	c1->cd(3);
	hupRes->Fit("gaus");
	hupRes->Draw();

	c1->cd(4);
	hvpRes->Fit("gaus");
	hvpRes->Draw();

	c1->cd(5);
	huRes->Fit("gaus");
	huRes->Draw();

	c1->cd(6);
	hvRes->Fit("gaus");
	hvRes->Draw();

	c1->Write();

	TCanvas* c2 = new TCanvas();
	c2->Divide(2, 3);

	c2->cd(1);
	hqopPu->Fit("gaus");
	hqopPu->Draw();

//	c2->cd(2);
//	pVal->Fit("pol1");
//	pVal->Draw();

//	c2->cd(2);
////		hchi2Forward->Fit("");
//	hchi2Forward->Draw();

	c2->cd(3);
	hupPu->Fit("gaus");
	hupPu->Draw();

	c2->cd(4);
	hvpPu->Fit("gaus");
	hvpPu->Draw();

	c2->cd(5);
	huPu->Fit("gaus");
	huPu->Draw();

	c2->cd(6);
	hvPu->Fit("gaus");
	hvPu->Draw();

	c2->Write();
	delete fitter;

	// open event display
	display->open();

}

