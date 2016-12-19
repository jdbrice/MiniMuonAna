#ifndef ST_MINI_MC_ANALYZER_H
#define ST_MINI_MC_ANALYZER_H

#include "TreeAnalyzer.h"
#include "CutCollection.h"
using namespace jdb;

#include "TLorentzVector.h"

#include "StMiniMcReader.h"
#include <memory>

class StMiniMcAnalyzer : public TreeAnalyzer
{
public:
	virtual const char* classname() const { return "StMiniMcAnalyzer"; }
	StMiniMcAnalyzer(){}
	~StMiniMcAnalyzer(){}

	virtual void initialize(){
		TreeAnalyzer::initialize();

		rdr = shared_ptr<StMiniMcReader>( new StMiniMcReader( chain ) );
		book->cd();

		mcTrackCuts.init( config, nodePath + ".McTrackCuts" );
		mcTrackCuts.setDefault( "pT", 0.0, 1000 );
		mcTrackCuts.setDefault( "eta", -10.0, 10.0 );
		mcTrackCuts.setDefault( "parentGeantId", -1, 100000 );

		mcTrackCuts.report();

	}

protected:
	shared_ptr<StMiniMcReader> rdr;
	int nPosMuon = 0;
	int nNegMuon = 0;
	vector<StTinyMcTrack*> posMuons;
	vector<StTinyMcTrack*> negMuons;
	CutCollection mcTrackCuts;

	virtual bool keepEvent(){
		return true;
	}
	virtual void analyzeEvent(){


		book->fill( "mCentrality", rdr->event->centrality() );
		book->fill( "mCentralMult", rdr->event->centralMult() );


		nPosMuon = 0;
		nNegMuon = 0;

		TClonesArray * mcTracks, *rcTracks;
		mcTracks = rdr->event->tracks( (Category)MC );
		rcTracks = rdr->event->tracks( (Category)MATCHED );

		TIter iMC(mcTracks);
		StTinyMcTrack * mcTrack;
		
		posMuons.clear();
		negMuons.clear();

		while( (mcTrack=(StTinyMcTrack*) iMC.Next()) ) {
			if ( keepMcTrack( mcTrack ) )
				analyzeMcTrack( mcTrack );
		} // while

		// INFO( classname(), "pos=" << nPosMuon << ", neg=" << nNegMuon );
		book->fill( "mc_n_pos_vs_neg", nPosMuon, nNegMuon );
		if ( nPosMuon > 0 && nNegMuon > 0 )
			book->fill( "mu_charge_ratio", nPosMuon / (float)nNegMuon );


		for ( int iPos = 0; iPos < posMuons.size(); iPos++ ){
			// unlike sign pairs
			for ( int iNeg = 0; iNeg < negMuons.size(); iNeg++ ){
				analyzeMcPair( posMuons[ iPos ], negMuons[ iNeg ] );
				if ( iPos == 0 ){
					for ( int iNeg2 = iNeg; iNeg2 < negMuons.size(); iNeg2++ ){
						if ( iNeg == iNeg2 ) continue;
						analyzeMcPair( negMuons[ iNeg ], negMuons[ iNeg2 ] );	
					}
				}
			}	

			for ( int iPos2 = iPos; iPos2 < posMuons.size(); iPos2++ ){
				if ( iPos == iPos2 ) continue;
				analyzeMcPair( posMuons[ iPos ], posMuons[ iPos2 ] );
			}
		}




		TIter iRC(rcTracks);
		StMiniMcPair * rcTrack;

		while( (rcTrack=(StMiniMcPair*) iRC.Next()) ) {
			analyzeRcTrack( rcTrack );
		} // while
	}

	virtual bool keepMcTrack( StTinyMcTrack* _track ){

		if (_track->ptMc() < mcTrackCuts[ "pT" ]->min )
			return false;

		if ( !mcTrackCuts["eta"]->inInclusiveRange( _track->etaMc() ) )
			return false;

		if ( !mcTrackCuts["parentGeantId"]->inInclusiveRange( _track->parentGeantId() ) )
			return false;

		return true;
	}


	virtual void analyzeMcPair( StTinyMcTrack* _trackA, StTinyMcTrack* _trackB ){
		
		if ( _trackA->key() == _trackB->key() ) return;
		int gIdA = _trackA->geantId();
		int pgIdA = _trackA->parentGeantId();
		int gIdB = _trackB->geantId();
		int pgIdB = _trackB->parentGeantId();
		
		TLorentzVector lv1, lv2, lv;
		lv1.SetXYZM( _trackA->pxMc(), _trackA->pyMc(), _trackA->pzMc(), 0.1057 );
		lv2.SetXYZM( _trackB->pxMc(), _trackB->pyMc(), _trackB->pzMc(), 0.1057 );
		lv = lv1 + lv2;

		if ( (gIdA == 5 && gIdB == 6) || (gIdA == 6 && gIdB == 5) ){
			book->fill( "mc_uls_mass", lv.M() );
			book->fill( "mc_uls_mass_vs_pt", lv.M(), lv.Pt() );

			if ( _trackA->ptMc() > 0.5 && _trackB->ptMc() > 0.5 )
				book->fill( "mc_uls_mass_0p5", lv.M() );
			if ( _trackA->ptMc() > 0.6 && _trackB->ptMc() > 0.6 )
				book->fill( "mc_uls_mass_0p6", lv.M() );
			if ( _trackA->ptMc() > 0.7 && _trackB->ptMc() > 0.7 )
				book->fill( "mc_uls_mass_0p7", lv.M() );
			if ( _trackA->ptMc() > 0.8 && _trackB->ptMc() > 0.8 )
				book->fill( "mc_uls_mass_0p8", lv.M() );
			if ( _trackA->ptMc() > 0.9 && _trackB->ptMc() > 0.9 )
				book->fill( "mc_uls_mass_0p9", lv.M() );
			if ( _trackA->ptMc() > 1.0 && _trackB->ptMc() > 1.0 )
				book->fill( "mc_uls_mass_1p0", lv.M() );



			if ( _trackA->parentKey() == _trackB->parentKey() ){
				book->fill( "mc_signal_uls_mass", lv.M() );
				book->fill( "mc_signal_uls_parents", pgIdA );
			}
		}

		if ( gIdA == 5 && gIdB == 5 ){

			// if ( lv.M()  < 0.215 ){
			// 	INFO( classname(), "trackKeyA: " << _trackA->key() << " trackKeyB: " << _trackB->key() )
			// 	INFO( classname(), "parentA : " << pgIdA << "parentKeyA: " << _trackA->parentKey());
			// 	INFO( classname(), "parentB : " << pgIdB << "parentKeyB: " << _trackB->parentKey() );
			// }
			book->fill( "mc_ls_mass", lv.M() );	
			book->fill( "mc_ls_pos_mass", lv.M() );

			book->fill( "mc_ls_mass_vs_pt", lv.M(), lv.Pt() );	
			book->fill( "mc_ls_pos_mass_vs_pt", lv.M(), lv.Pt() );

			if ( _trackA->parentKey() == _trackB->parentKey() ){
				book->fill( "mc_signal_ls_mass", lv.M() );
				book->fill( "mc_signal_ls_parents", pgIdA );
			}
		}
		if ( gIdA == 6 && gIdB == 6 ){
			book->fill( "mc_ls_mass", lv.M() );	
			book->fill( "mc_ls_neg_mass", lv.M() );

			book->fill( "mc_ls_mass_vs_pt", lv.M(), lv.Pt() );	
			book->fill( "mc_ls_neg_mass_vs_pt", lv.M(),lv.Pt() );

			if ( _trackA->parentKey() == _trackB->parentKey() ){
				book->fill( "mc_signal_ls_mass", lv.M() );
				book->fill( "mc_signal_ls_parents", pgIdA );
			}
		}


	}

	virtual void analyzeMcTrack( StTinyMcTrack* _track ) {
		int gId = _track->geantId();
		int pgId = _track->parentGeantId();
		
		book->fill( "mc_plcs", gId );
		book->fill( "mc_parents", pgId );

		if ( gId == 5 ){
			book->fill( "mup_mc_pT", _track->ptMc() );
			book->fill( "mup_mc_parents", pgId );
			posMuons.push_back( _track );
			nPosMuon++;

		} else if ( gId == 6 ){
			book->fill( "mum_mc_pT", _track->ptMc() );
			book->fill( "mum_mc_parents", pgId );
			negMuons.push_back( _track );
			nNegMuon ++;
		}
		
	}

	virtual void analyzeRcTrack( StMiniMcPair *_track ){
		book->fill( "rc_pT", _track->ptPr() );
	}


	void addGEANTLabels( TH1* h ){
		TAxis * x = h->GetXaxis();


		x->SetBinLabel( 1, "No ID" );
	    x->SetBinLabel( 2, "#gamma");
	    x->SetBinLabel( 3, "e^{+}" );
	    x->SetBinLabel( 4, "e^{-}" );
	    x->SetBinLabel( 5, "#nu" );
	    x->SetBinLabel( 6, "#mu^{+}" );
	    x->SetBinLabel( 7, "#mu^{-}" );
	    x->SetBinLabel( 8, "#pi" );
	    x->SetBinLabel( 9, "#pi^{+}" );
	    x->SetBinLabel( 10, "#pi^{-}" );
	    x->SetBinLabel( 11, "K^{0}_{L}" );
	    x->SetBinLabel( 12, "K^{+}" );
	    x->SetBinLabel( 13, "K^{-}" );
	    x->SetBinLabel( 14, "n" );
	    x->SetBinLabel( 15, "P" );
	    x->SetBinLabel( 16, "#bar{P}" );
	    x->SetBinLabel( 17, "K^{0}_{S}" );
	    x->SetBinLabel( 18, "#eta" );
	    x->SetBinLabel( 19, "#Lambda" );
	    x->SetBinLabel( 20, "#Sigma^{+}" );
	    x->SetBinLabel( 21, "#Sigma^{0}" );
	    x->SetBinLabel( 22, "#Sigma^{-}" );
	    x->SetBinLabel( 23, "#Xi^{0}" );
	    x->SetBinLabel( 24, "#Xi^{-}" );
	    x->SetBinLabel( 25, "#Omega^{-}" );
	    x->SetBinLabel( 26, "#bar{n}" );
	    x->SetBinLabel( 27, "#bar{#Lambda}" );
	    x->SetBinLabel( 28, "#bar{#Sigma^{-}}" );
	    x->SetBinLabel( 29, "#bar{#Sigma^{0}}" );
	    x->SetBinLabel( 30, "#bar{#Sigma^{+}}" );
	    x->SetBinLabel( 31, "#bar{#Xi^{0}}" );
	    x->SetBinLabel( 32, "#bar{#Xi^{+}}" );
	    x->SetBinLabel( 33, "#bar{#Omega^{+}}" );
	}

	virtual void postEventLoop(){
		addGEANTLabels( book->get( "mc_plcs" ) );
		addGEANTLabels( book->get( "mc_parents" ) );

		addGEANTLabels( book->get( "mup_mc_parents" ) );
		addGEANTLabels( book->get( "mum_mc_parents" ) );
	}


	
};



#endif