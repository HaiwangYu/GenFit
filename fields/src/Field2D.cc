/* Copyright 2008-2010, Technische Universitaet Muenchen,
   Authors: Christian Hoeppner & Sebastian Neubert & Johannes Rauch

   This file is part of GENFIT.

   GENFIT is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published
   by the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   GENFIT is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with GENFIT.  If not, see <http://www.gnu.org/licenses/>.
*/
#include "Field2D.h"

#include <TFile.h>
#include <TNtuple.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TH2D.h>

#include <iostream>
#include <cmath>

#define LogDEBUG std::cout<<"DEBUG: "<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<"\n"

namespace genfit {

bool Field2D::initialize(std::string inname) {
	TFile *fin = TFile::Open(inname.data(),"READ");
	if(!fin)
	{
		//TODO Output some msg
		LogDEBUG;
		return false;
	}
	TNtuple *T = (TNtuple *) fin->Get("fieldmap");
	if(!T)
	{
		//TODO Output some msg
		LogDEBUG;
		return false;
	}

	Float_t r;
	Float_t z;
	Float_t br;
	Float_t bz;
	T->SetBranchAddress("r",&r);
	T->SetBranchAddress("z",&z);
	T->SetBranchAddress("br",&br);
	T->SetBranchAddress("bz",&bz);

	for(long ientry=0;ientry<T->GetEntries();ientry++){
		//if (T->GetEntry(ientry) < 0) continue;
		T->GetEntry(ientry);
		field_map_.insert(BPAIR2D(boost::make_tuple(r,z),boost::make_tuple(br*10.,bz*10.))); // Tesla => kGauss
	}
	fin->Close();

	return true;
}

void Field2D::plot(std::string option){

	std::cout<<"field_map_.size():"<<field_map_.size()<<"\n";

	TH2D *hbr = new TH2D("hbr","|B_{r}| [kGauss]; z [cm]; r [cm]",800, -400,400,300, 0, 300);
	TH2D *hbz = new TH2D("hbz","|B_{z}| [kGauss]; z [cm]; r [cm]",800, -400,400,300, 0, 300);
	for(BMAP2D::iterator it = field_map_.begin(); it != field_map_.end(); it++)
	{
		double r = it->first.get<0>();
		double z = it->first.get<1>();
		double br = it->second.get<0>();
		double bz = it->second.get<1>();
//
//		std::cout<<
//				r <<" ,"<<
//				z <<" ,"<<
//				br<<" ,"<<
//				bz<<"\n";

		hbr->Fill(z,r,abs(br));
		hbz->Fill(z,r,abs(bz));
	}

	TCanvas *c0 = new TCanvas("c0","c0");
	c0->Divide(1,2);
	c0->cd(1);
	hbr->SetStats(0);
	hbr->Draw("colz");
	c0->cd(2);
	hbz->SetStats(0);
	hbz->Draw("colz");

	c0->SaveAs("Field2D_plot.root");
	c0->SaveAs("Field2D_plot.png");

	delete c0;
	delete hbr;
	delete hbz;


	return;
}

TVector3 Field2D::get(const TVector3& v) const {
  double x = v.x();
  double y = v.y();
  double z = v.z();
  double Bx;
  double By;
  double Bz;
  get(x,y,z,Bx,By,Bz);
  return TVector3(Bx,By,Bz);
}

void Field2D::get(const double&x, const double&y, const double&z, double& Bx, double& By, double& Bz) const {
	double r = sqrt(x*x + y*y);
	boost::tuple<double, double> Brz = field_map_.find(boost::make_tuple(r,z))->second;

	Bx = x/r*Brz.get<0>();
	By = y/r*Brz.get<0>();
	Bz = Brz.get<1>();
}

} /* End of namespace genfit */
