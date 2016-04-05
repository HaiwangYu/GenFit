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
//#include <boost/shared_ptr.hpp>
//#include <boost/make_shared.hpp>
#include <KalmanFittedStateOnPlane.h>
//#include <AbsKalmanFitter.h>
//#include <KalmanFitter.h>
//#include <KalmanFitterRefTrack.h>
#include <KalmanFitterInfo.h>
//#include <KalmanFitStatus.h>
#include <KalmanFitter.h>


#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TProfile.h>
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
	//double resolution_detector_xy = 0.1;
	double resolution_detector_xy = 0.005/3.; //50 micron
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
	fieldMap->re_scale(1.4/1.5);// Re-scale to 1.4 T
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
	TTree *Tout = new TTree("T","Test GenFit2");
	//fout->cd();

	TH1D *hmomRes = new TH1D("hmomRes", "mom residual", 500,
			-200 * resolution_detector_xy * momentum_initial / nMeasurements,
			200 * resolution_detector_xy * momentum_initial / nMeasurements);
	TH1D *hupRes = new TH1D("hupRes", "u' residual", 500,
			-15 * resolution_detector_xy / nMeasurements,
			15 * resolution_detector_xy / nMeasurements);
	TH1D *hvpRes = new TH1D("hvpRes", "v' residual", 500,
			-15 * resolution_detector_xy / nMeasurements,
			15 * resolution_detector_xy / nMeasurements);
	TH1D *huRes = new TH1D("huRes", "u residual", 500, -15 * resolution_detector_xy,
			15 * resolution_detector_xy);
	TH1D *hvRes = new TH1D("hvRes", "v residual", 500, -15 * resolution_detector_xy,
			15 * resolution_detector_xy);

	TH1D *hqopPu = new TH1D("hqopPu", "q/p pull", 200, -6., 6.);
//	TH1D *pVal = new TH1D("pVal", "p-value", 100, 0., 1.00000001);
//	pVal->SetMinimum(0);
	TH1D *hupPu = new TH1D("hupPu", "u' pull", 200, -6., 6.);
	TH1D *hvpPu = new TH1D("hvpPu", "v' pull", 200, -6., 6.);
	TH1D *huPu = new TH1D("huPu", "u pull", 200, -6., 6.);
	TH1D *hvPu = new TH1D("hvPu", "v pull", 200, -6., 6.);

	TH1D *hchi2Forward = new TH1D("hchi2Forward", "Fit #chi^{2}", 100, 0, 3.);

	//TH2D *hpT_residual_vs_pT = new TH2D("hpT_residual_vs_pT", "#delta pT/pT; pT[GeV/c]; #Delta pT/pT", 200, 1, 40, 100, 0, 0.1);
	TH2D *hpT_residual_vs_pT = new TH2D("hpT_residual_vs_pT", "#Delta pT/pT; pT[GeV/c]; #Delta pT/pT", 80, 0.5, 40.5, 1000, -1, 1);
	TH2D *hADpT_residual_vs_pT = (TH2D*) hpT_residual_vs_pT->Clone("hADpT_residual_vs_pT");


	TH1D *hpx_residual = new TH1D("hpx_residual", "px^{GenFit} - px^{True}; #Delta px [GeV/c]; ", 100, -0.5, 0.5);
	TH1D *hpy_residual = new TH1D("hpy_residual", "py^{GenFit} - py^{True}; #Delta py [GeV/c]; ", 100, -0.5, 0.5);
	TH1D *hpz_residual = new TH1D("hpz_residual", "pz^{GenFit} - pz^{True}; #Delta pz [GeV/c]; ", 100, -0.5, 0.5);
	TH1D *hptot_residual = new TH1D("hptot_residual", "ptot^{GenFit} - ptot^{True}; #Delta ptot [GeV/c]; ", 100, -1, 1);

	TH1D *hpT_residual = new TH1D("hpT_residual", "p_{T}^{GenFit} - p_{T}^{True}; #Delta p{T} [GeV/c]; ", 100, -0.5, 0.5);
	TH1D *hpT_diff = new TH1D("hpT_diff", "p_{T}^{GenFit} - p_{T}^{Current}; #Delta p{T} [GeV/c]; ", 100, -0.5, 0.5);

	TH2D *hpT_resolution_vs_diff = new TH2D("hpT_resolution_vs_diff","", 100, -0.5, 0.5, 100, -0.5, 0.5);



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


	Float_t Cluster_x[NLAYERS];
	Float_t Cluster_y[NLAYERS];
	Float_t Cluster_z[NLAYERS];
	Float_t Cluster_size_dphi[NLAYERS];
	Float_t Cluster_size_dz[NLAYERS];
	Float_t True_px;
	Float_t True_py;
	Float_t True_pz;
	Float_t AlanDion_px;
	Float_t AlanDion_py;
	Float_t AlanDion_pz;
	Float_t AlanDion_dca2d;
	Int_t nhits;

	T->SetBranchAddress("nhits", &nhits);
	T->SetBranchAddress("gpx", &True_px);
	T->SetBranchAddress("gpy", &True_py);
	T->SetBranchAddress("gpz", &True_pz);
	T->SetBranchAddress("px", &AlanDion_px);
	T->SetBranchAddress("py", &AlanDion_py);
	T->SetBranchAddress("pz", &AlanDion_pz);
	T->SetBranchAddress("dca2d", &AlanDion_dca2d);
	T->SetBranchAddress("x", Cluster_x);
	T->SetBranchAddress("y", Cluster_y);
	T->SetBranchAddress("z", Cluster_z);
	T->SetBranchAddress("size_dphi", Cluster_size_dphi);
	T->SetBranchAddress("size_dz", Cluster_size_dz);

	Float_t GenFit_px;
	Float_t GenFit_py;
	Float_t GenFit_pz;
	Float_t GenFit_dca2d;
	Float_t GenFit_chi2_ndf;

	Tout->Branch("GenFit_px",&GenFit_px,"GenFit_px/F");
	Tout->Branch("GenFit_py",&GenFit_py,"GenFit_py/F");
	Tout->Branch("GenFit_pz",&GenFit_pz,"GenFit_pz/F");
	Tout->Branch("GenFit_dca2d",&GenFit_dca2d,"GenFit_dca2d/F");
	Tout->Branch("GenFit_chi2_ndf",&GenFit_chi2_ndf,"GenFit_chi2_ndf/F");

	Tout->Branch("AlanDion_px",&AlanDion_px,"AlanDion_px/F");
	Tout->Branch("AlanDion_py",&AlanDion_py,"AlanDion_py/F");
	Tout->Branch("AlanDion_pz",&AlanDion_pz,"AlanDion_pz/F");
	Tout->Branch("AlanDion_dca2d",&AlanDion_dca2d,"AlanDion_dca2d/F");
	Tout->Branch("True_px",&True_px,"True_px/F");
	Tout->Branch("True_py",&True_py,"True_py/F");
	Tout->Branch("True_pz",&True_pz,"True_pz/F");


	// main loop
	for (unsigned int ientry = 0; ientry < T->GetEntries(); ++ientry) {
	//for (unsigned int ientry = 0; ientry < 100000; ++ientry) {
		//T->GetEntry(atoi(argv[1]));
		T->GetEntry(ientry);

		if(nhits < 0)
		{
			LogDEBUG;
			continue;
		}

		// true start values
		TVector3 init_pos(0, 0, 0); //cm
		TVector3 True_mom(True_px, True_py, True_pz);

		const int pdg = -13; //-13: mu+, 13: mu-
		const double charge =
				TDatabasePDG::Instance()->GetParticle(pdg)->Charge() / (3.);

		// Seed: use smeared values
		const bool smearPosMom = true; // init the Reps with smeared init_pos and True_mom
		const double posSmear = 10 * resolution_detector_xy;     // cm
		const double momSmear = 3. / 180. * TMath::Pi();     // rad
		const double momMagSmear = 0.1;   // relative

		TVector3 seed_pos(init_pos);
		TVector3 seed_mom(True_mom);
		if (smearPosMom) {
			seed_pos.SetX(gRandom->Gaus(seed_pos.X(), posSmear));
			seed_pos.SetY(gRandom->Gaus(seed_pos.Y(), posSmear));
			seed_pos.SetZ(gRandom->Gaus(seed_pos.Z(), posSmear));

			seed_mom.SetPhi(gRandom->Gaus(True_mom.Phi(), momSmear));
			seed_mom.SetTheta(gRandom->Gaus(True_mom.Theta(), momSmear));
			seed_mom.SetMag(
					gRandom->Gaus(True_mom.Mag(),
							momMagSmear * True_mom.Mag()));
		}

		// approximate covariance
		TMatrixDSym seed_cov(6);

//		measurementCreator.setResolution(resolution_detector_xy);

		for (int ilayer = 0; ilayer < 3; ++ilayer)
			seed_cov(ilayer, ilayer) = resolution_detector_xy
					* resolution_detector_xy;
		for (int ilayer = 3; ilayer < 6; ++ilayer)
			seed_cov(ilayer, ilayer) = pow(
					resolution_detector_xy / nMeasurements / sqrt(3), 2);

		// trackrep
		genfit::AbsTrackRep* rep = new genfit::RKTrackRep(pdg);
		genfit::AbsTrackRep* measurementCreatorRep = new genfit::RKTrackRep(
				pdg);

		// seed start state
		genfit::MeasuredStateOnPlane seedMSoP(rep);
		seedMSoP.setPosMomCov(seed_pos, seed_mom, seed_cov);
		const genfit::StateOnPlane seedSoP(seedMSoP);

		genfit::MeasuredStateOnPlane initMSoP(measurementCreatorRep);
		initMSoP.setPosMomCov(init_pos, True_mom, seed_cov);
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

//		genfit::SharedPlanePtr destPlane = boost::make_shared<genfit::DetPlane>(
//				TVector3(0, 0, dest_z), TVector3(0, 0, 1));
		genfit::SharedPlanePtr destPlane(new genfit::DetPlane(
				TVector3(0, 0, dest_z), TVector3(0, 0, 1)));

		int measurementCounter_ = 0;
		// create smeared measurements and add to track
		try {
			for (unsigned int ilayer = 0; ilayer < NLAYERS; ++ilayer) {

				genfit::SharedPlanePtr plane(
						new genfit::DetPlane(
								TVector3(Cluster_x[ilayer], Cluster_y[ilayer],
										Cluster_z[ilayer]),
								TVector3(Cluster_x[ilayer], Cluster_y[ilayer],
										0)));

//				std::cout << "DEBUG: " << " position (cm):" <<
//						Cluster_x[ilayer] << ","<<
//						Cluster_y[ilayer] << ","<<
//						Cluster_z[ilayer] << "\n";

				int nDim = 2;
				TVectorD hitCoords(nDim);
				TMatrixDSym hitCov(nDim);

				hitCoords(0) = 0;
				hitCoords(1) = 0;

				//double resolution_detector_z = 0.0425/3;
				hitCov(0, 0) = Cluster_size_dphi[ilayer]*Cluster_size_dphi[ilayer]/12.;
				hitCov(1, 1) = Cluster_size_dz[ilayer]*Cluster_size_dz[ilayer]/12.;

//				LogDEBUG;
//				plane->Print();

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
		//fitter->setDebugLvl(1);
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
			if((fabs(AlanDion_px) < 100))
				std::cout<<"IMPORTANT: Not reco'ed by GenFit: "<<ientry<<"\n";
			continue;
		}


//		fitTrack.Print();
		double chi2 = fitTrack.getFitStatus(rep)->getChi2();
		double ndf = fitTrack.getFitStatus(rep)->getNdf();
		hchi2Forward->Fill(chi2 / ndf);

//! Naiive extrapolate, which seems to be wrong.
//		fitTrack.getCardinalRep()->extrapolateToPoint(currentState,
//				TVector3(0, 0, 0));
//		std::cout << "DEBUG: extrapolateToPoint(0,0,0)\n";
//		fitTrack.getCardinalRep()
//		rep->extrapolateToPlane(currentState, seedSoP.getPlane());
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
//			rep->extrapolateToPlane(kfsop, initSoP.getPlane());
			rep->extrapolateToLine(kfsop, TVector3(0.,0.,0.), TVector3(0.,0.,1.));
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

		hmomRes->Fill((charge / state[0] - True_mom.Mag())/True_mom.Mag());
		hupRes->Fill((state[1] - referenceState[1]));
		hvpRes->Fill((state[2] - referenceState[2]));
		huRes->Fill((state[3] - referenceState[3]));
		hvRes->Fill((state[4] - referenceState[4]));

		hqopPu->Fill((state[0] - referenceState[0]) / sqrt(cov[0][0]));
		hupPu->Fill((state[1] - referenceState[1]) / sqrt(cov[1][1]));
		hvpPu->Fill((state[2] - referenceState[2]) / sqrt(cov[2][2]));
		huPu->Fill((state[3] - referenceState[3]) / sqrt(cov[3][3]));
		hvPu->Fill((state[4] - referenceState[4]) / sqrt(cov[4][4]));


		//kfsop.Print();

		TVector3 GenFit_mom = kfsop.getMom();

		TVector3 AlanDion_mom(AlanDion_px,AlanDion_py,AlanDion_pz);

		hpT_residual_vs_pT->Fill(True_mom.Pt(),(GenFit_mom.Pt() - True_mom.Pt())/True_mom.Pt());
		hADpT_residual_vs_pT->Fill(True_mom.Pt(),(AlanDion_mom.Pt() - True_mom.Pt())/True_mom.Pt());


		hpx_residual->Fill(GenFit_mom.Px() - True_mom.Px());
		hpy_residual->Fill(GenFit_mom.Py() - True_mom.Py());
		hpz_residual->Fill(GenFit_mom.Pz() - True_mom.Pz());

		hpT_residual->Fill(GenFit_mom.Pt() - True_mom.Pt());
		hptot_residual->Fill(GenFit_mom.Mag() - True_mom.Mag());

		hpT_diff->Fill(GenFit_mom.Pt()- AlanDion_mom.Pt());

		if(!(fabs(AlanDion_px) < 100))
			std::cout<<"IMPORTANT: Not reco'd by AD: "<<ientry<<"\n";

		hpT_resolution_vs_diff->Fill(GenFit_mom.Pt()- AlanDion_mom.Pt(),GenFit_mom.Mag() - True_mom.Mag());


		GenFit_px = GenFit_mom.Px();
		GenFit_py = GenFit_mom.Py();
		GenFit_pz = GenFit_mom.Pz();
		GenFit_chi2_ndf = chi2 / ndf;
		GenFit_dca2d = state[3]; //u
		Tout->Fill();
	}    // end loop over events

	fout->cd();
	Tout->Write();


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

	c2->cd(2);
//hchi2Forward->Fit("");
	hchi2Forward->Draw();

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


	TF1 *tf_pT_resolution = new TF1("tf_pT_resolution","sqrt([0]*[0] + x*x*[1]*[1])", 0, 40);
	tf_pT_resolution->SetParameters(0,0);
	TCanvas *c3 = new TCanvas("c3","c3");
	c3->Divide(2,1);
	c3->cd(1);
//	hpT_residual_vs_pT->Draw("colz");
//	c3->cd(2);
	hpT_residual_vs_pT->FitSlicesY();
	TH1D *hpT_resolution_vs_pT = (TH1D*)gDirectory->Get("hpT_residual_vs_pT_2");
	hpT_resolution_vs_pT->SetTitle("GenFit: #sigma_{p_{T}}/p_{T}; p_{T}[GeV/c]; #sigma_{p_{T}}/p_{T}");
	hpT_resolution_vs_pT->SetMarkerStyle(20);
	hpT_resolution_vs_pT->Draw("e");
	hpT_resolution_vs_pT->Fit(tf_pT_resolution);
	c3->cd(2);
//	hADpT_residual_vs_pT->Draw("colz");
//	c3->cd(4);
	hADpT_residual_vs_pT->FitSlicesY();
	TH1D *hADpT_resolution_vs_pT = (TH1D*)gDirectory->Get("hADpT_residual_vs_pT_2");
	hADpT_resolution_vs_pT->SetTitle("Current Kalman: #sigma_{p_{T}}/p_{T}; p_{T}[GeV/c]; #sigma_{p_{T}}/p_{T}");
	hADpT_resolution_vs_pT->SetMarkerStyle(20);
	hADpT_resolution_vs_pT->Draw("e");
	hADpT_resolution_vs_pT->Fit(tf_pT_resolution);
	c3->Update();

	TCanvas *c4 = new TCanvas("c4","c4");
	c4->Divide(2,2);
	c4->cd(1);
	//hpx_residual->Draw();
	hpT_diff->Draw();
	c4->cd(2);
	hpT_residual->Draw();
	c4->cd(3);
	//hpz_residual->Draw();
	hpT_resolution_vs_diff->SetTitle(";p_{T}^{GenFit} - p_{T}^{Current} [GeV/c];p_{T}^{GenFit} - p_{T}^{True} [GeV/c]");
	hpT_resolution_vs_diff->Draw("colz");
	hpT_resolution_vs_diff->ProfileX()->Draw("same");
	c4->cd(4);
	hptot_residual->Draw();

	c4->Update();

	fout->cd();
	c1->Write();
	c2->Write();
	c3->Write();
	c4->Write();

	fout->Close();
	fPHG4Hits->Close();


	delete fitter;

	// open event display
	display->open();

}

