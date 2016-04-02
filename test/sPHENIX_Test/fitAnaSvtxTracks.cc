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
#include <KalmanFitter.h>



#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <TStyle.h>

#include <Field2D.h>

#define LogDEBUG    std::cout<<"DEBUG: "<<__LINE__<<"\n"


#define NLAYERS 7

int main(int argc, char **argv) {

	double magnetic_field = 20.; //kGauss
	double momentum_initial = 10.; // GeV

	double dest_z = 10; //cm
	std::cout << magnetic_field << "\n";
	//double resolution_detector = 0.1;
	double resolution_detector = 0.005; //50 micron
	unsigned int nMeasurements = NLAYERS; //
	gRandom->SetSeed(14);

	// init MeasurementCreator
	//genfit::MeasurementCreator measurementCreator;

	// init geometry and mag. field
//	new TGeoManager("Geometry", "Geane geometry");
//	TGeoManager::Import("genfitGeom.root");
	new TGeoManager("Default", "Geane geometry");
	TGeoManager::Import("sPHENIX_Geo.root");
//	genfit::FieldManager::getInstance()->init(
//			new genfit::ConstField(0., 0., magnetic_field)); // kGauss
	genfit::Field2D *fieldMap = new genfit::Field2D("sPHENIX.2d.root");
	double bx,by,bz;
	fieldMap->get(1,0,0,bx,by,bz);
	std::cout<<"DEBUG: "<<bx<<","<<by<<","<<bz<<"\n";
	//fieldMap->plot();
	genfit::FieldManager::getInstance()->init(
			fieldMap);
	genfit::MaterialEffects::getInstance()->init(
			new genfit::TGeoMaterialInterface());

	// init event display
	genfit::EventDisplay* display = genfit::EventDisplay::getInstance();

	// init fitter
	genfit::AbsKalmanFitter* fitter = new genfit::KalmanFitterRefTrack();

	TFile *fout = TFile::Open("results.root","recreate");
	fout->cd();

	TH1D *hmomRes = new TH1D("hmomRes", "mom residual", 500,
			-200 * resolution_detector * momentum_initial / nMeasurements,
			200 * resolution_detector * momentum_initial / nMeasurements);
	TH1D *hupRes = new TH1D("hupRes", "u' residual", 500,
			-15 * resolution_detector / nMeasurements,
			15 * resolution_detector / nMeasurements);
	TH1D *hvpRes = new TH1D("hvpRes", "v' residual", 500,
			-15 * resolution_detector / nMeasurements,
			15 * resolution_detector / nMeasurements);
	TH1D *huRes = new TH1D("huRes", "u residual", 500, -15 * resolution_detector,
			15 * resolution_detector);
	TH1D *hvRes = new TH1D("hvRes", "v residual", 500, -15 * resolution_detector,
			15 * resolution_detector);

	TH1D *hqopPu = new TH1D("hqopPu", "q/p pull", 200, -6., 6.);
//	TH1D *pVal = new TH1D("pVal", "p-value", 100, 0., 1.00000001);
//	pVal->SetMinimum(0);
	TH1D *hupPu = new TH1D("hupPu", "u' pull", 200, -6., 6.);
	TH1D *hvpPu = new TH1D("hvpPu", "v' pull", 200, -6., 6.);
	TH1D *huPu = new TH1D("huPu", "u pull", 200, -6., 6.);
	TH1D *hvPu = new TH1D("hvPu", "v pull", 200, -6., 6.);

	TH1D *hchi2Forward = new TH1D("hchi2Forward", "Fit #chi^{2}", 100, 0, 3.);

	TH2D *hmomresolution = new TH2D("hmomresolution", "#delta p/p; 1/p[1/(GeV/c)]; #delta p/p", 200, 0.05, 2, 100, 0, 0.1);
	TH2D *hpT_resolution = new TH2D("hpT_resolution", "#delta p_{T}/p_{T}; 1/p_{T} [1/(GeV/c)]; #delta p_{T}/p_{T}", 200, 0.05, 2, 100, 0, 0.1);
	TH1D *hpx_residual = new TH1D("hpx_residual", "px^{GenFit} - px^{Current}; #Delta px [GeV/c]; ", 100, -2, 2);
	TH1D *hpy_residual = new TH1D("hpy_residual", "py^{GenFit} - py^{Current}; #Delta py [GeV/c]; ", 100, -2, 2);
	TH1D *hpz_residual = new TH1D("hpz_residual", "pz^{GenFit} - pz^{Current}; #Delta pz [GeV/c]; ", 100, -2, 2);
	TH1D *hptot_residual = new TH1D("hptot_residual", "ptot^{GenFit} - ptot^{Current}; #Delta ptot [GeV/c]; ", 100, -2, 2);


	TFile *fPHG4Hits = TFile::Open("AnaSvtxTracksForGenFit.root", "read");
	if (!fPHG4Hits) {
		std::cout << "No TFile Openned: " << __LINE__ << "\n";
		return -1;
	}
	TTree *T = (TTree*) fPHG4Hits->Get("tracks");
	if (!T) {
		std::cout << "No TTree Found: " << __LINE__ << "\n";
		return -1;
	}


	Float_t reco_x[NLAYERS];
	Float_t reco_y[NLAYERS];
	Float_t reco_z[NLAYERS];
	Float_t true_px;
	Float_t true_py;
	Float_t true_pz;
	Int_t nhits;

	T->SetBranchAddress("nhits", &nhits);
	T->SetBranchAddress("px", &true_px);
	T->SetBranchAddress("py", &true_py);
	T->SetBranchAddress("pz", &true_pz);
	T->SetBranchAddress("x", reco_x);
	T->SetBranchAddress("y", reco_y);
	T->SetBranchAddress("z", reco_z);

	// main loop
	for (unsigned int ientry = 0; ientry < T->GetEntries(); ++ientry) {
	//for (unsigned int ientry = 0; ientry < 1; ++ientry) {
		//T->GetEntry(atoi(argv[1]));
		T->GetEntry(ientry);

		// true start values
		TVector3 init_pos(0, 0, 0); //cm
		TVector3 init_mom(true_px, true_py, true_pz);
		double true_p = init_mom.Mag();
		double true_pT = init_mom.Pt();

		const int pdg = -13; //-13: mu+, 13: mu-
		const double charge =
				TDatabasePDG::Instance()->GetParticle(pdg)->Charge() / (3.);

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

		for (int ilayer = 0; ilayer < 3; ++ilayer)
			seed_cov(ilayer, ilayer) = resolution_detector
					* resolution_detector;
		for (int ilayer = 3; ilayer < 6; ++ilayer)
			seed_cov(ilayer, ilayer) = pow(
					resolution_detector / nMeasurements / sqrt(3), 2);

		// trackrep
		genfit::AbsTrackRep* rep = new genfit::RKTrackRep(pdg);
		genfit::AbsTrackRep* measurementCreatorRep = new genfit::RKTrackRep(
				pdg);

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
		genfit::Track measurementCreatorTrack(measurementCreatorRep, initState,
				initCov);

		genfit::MeasuredStateOnPlane currentState = seedMSoP;

		genfit::SharedPlanePtr destPlane = boost::make_shared<genfit::DetPlane>(
				TVector3(0, 0, dest_z), TVector3(0, 0, 1));

		int measurementCounter_ = 0;
		// create smeared measurements and add to track
		try {
			for (unsigned int ilayer = 0; ilayer < NLAYERS; ++ilayer) {

				genfit::SharedPlanePtr plane(
						new genfit::DetPlane(
								TVector3(reco_x[ilayer], reco_y[ilayer],
										reco_z[ilayer]),
								TVector3(reco_x[ilayer], reco_y[ilayer],
										0)));

//				std::cout << "DEBUG: " << " position (cm):" <<
//						reco_x[ilayer] << ","<<
//						reco_y[ilayer] << ","<<
//						reco_z[ilayer] << "\n";

				int nDim = 2;
				TVectorD hitCoords(nDim);
				TMatrixDSym hitCov(nDim);

				hitCoords(0) = 0 + gRandom->Gaus(0, resolution_detector);
				hitCoords(1) = 0 + gRandom->Gaus(0, resolution_detector);

				for (int iDim = 0; iDim < nDim; iDim++) {
					hitCov(iDim, iDim) = resolution_detector
							* resolution_detector;
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
		fitter->setDebugLvl(1);
		fitter->processTrack(&fitTrack);

		//check
		assert(fitTrack.checkConsistency());

		if (ientry < 10) {
			// add track to event display
			display->addEvent(&fitTrack);
		}

		// check if fit was successful
		if (!fitTrack.getFitStatus(rep)->isFitConverged()) {
			std::cout
					<< "Track could not be fitted successfully! Fit is not converged! \n";
			//continue;
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
		//std::cout << "DEBUG: yuhw extrapolateToPlane: \n";
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
			//std::cout << "DEBUG: Official extrapolateToPlane: \n";
//			kfsop.Print();
		} catch (genfit::Exception& e) {
			std::cerr << "Exception, next track" << std::endl;
			std::cerr << e.what();
			continue;
		}

		TVectorD referenceState = initSoP.getState();
		TVectorD state = kfsop.getState();
		TMatrixDSym cov = kfsop.getCov();

		hmomRes->Fill((charge / state[0] - init_mom.Mag())/init_mom.Mag());
		//hmomRes->Fill((state[0] - 1./init_mom.Mag()));
		hupRes->Fill((state[1] - referenceState[1]));
		hvpRes->Fill((state[2] - referenceState[2]));
		huRes->Fill((state[3] - referenceState[3]));
		hvRes->Fill((state[4] - referenceState[4]));

		hqopPu->Fill((state[0] - referenceState[0]) / sqrt(cov[0][0]));
		hupPu->Fill((state[1] - referenceState[1]) / sqrt(cov[1][1]));
		hvpPu->Fill((state[2] - referenceState[2]) / sqrt(cov[2][2]));
		huPu->Fill((state[3] - referenceState[3]) / sqrt(cov[3][3]));
		hvPu->Fill((state[4] - referenceState[4]) / sqrt(cov[4][4]));

		hmomresolution->Fill(1./true_p,true_p*sqrt(cov[0][0]));


		//kfsop.Print();

		TVector3 reco_mom = kfsop.getMom();

//		double reco_p  = kfsop.getMom().Mag();
//		double reco_px = kfsop.getMom().Px();
//		double reco_py = kfsop.getMom().Py();
//		double reco_pz = kfsop.getMom().Pz();
//		double reco_pT = kfsop.getMom().Pt();

		hpx_residual->Fill(reco_mom.Px() - init_mom.Px());
		hpy_residual->Fill(reco_mom.Py() - init_mom.Py());
		hpz_residual->Fill(reco_mom.Pz() - init_mom.Pz());
		hptot_residual->Fill(reco_mom.Mag() - init_mom.Mag());


	}    // end loop over events
	fPHG4Hits->Close();


	gStyle->SetOptFit();
	fout->cd();

	TCanvas* c1 = new TCanvas("c1","c1");
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

	c1->Update();

	fout->cd();
	c1->Write();

	TCanvas* c2 = new TCanvas("c2","c2");
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

	c2->Update();
	fout->cd();
	c2->Write();

	TCanvas *c3 = new TCanvas("c3","c3");
	c3->cd();
	hmomresolution->Draw();
	c3->Update();

	TCanvas *c4 = new TCanvas("c4","c4");
	c4->Divide(2,2);
	c4->cd(1);
	hpx_residual->Draw();
	c4->cd(2);
	hpy_residual->Draw();
	c4->cd(3);
	hpz_residual->Draw();
	c4->cd(4);
	hptot_residual->Draw();

	c4->Update();


	delete fitter;

	// open event display
	display->open();

}

