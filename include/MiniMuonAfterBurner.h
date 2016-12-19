#ifndef MINI_MUON_AFTER_BURNER_H
#define MINI_MUON_AFTER_BURNER_H

#include "HistoAnalyzer.h"
using namespace jdb;

class MiniMuonAfterBurner : public HistoAnalyzer
{
public:
	virtual const char* classname() const { return "MiniMuonAfterBurner"; }
	MiniMuonAfterBurner(){}
	~MiniMuonAfterBurner(){}

	virtual void initialize(){
		HistoAnalyzer::initialize();
	}

protected:

	virtual void make(){

		TH1D * mc_ls_pos_mass = getH1D( "mc_ls_pos_mass" );
		TH1D * mc_ls_neg_mass = getH1D( "mc_ls_neg_mass" );
		
		book->cd();
		book->addClone( "pos_over_neg_vs_mass", mc_ls_pos_mass );

		book->get( "pos_over_neg_vs_mass" )->Divide( mc_ls_neg_mass );

	}
	
};

#endif