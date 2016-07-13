#include <algorithm>
#include <math.h>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/HeavyIonEvent/interface/EvtPlane.h"
#include "RecoHI/HiEvtPlaneAlgos/interface/HiEvtPlaneList.h"

#include "FWCore/Utilities/interface/RandomNumberGenerator.h"
#include "CLHEP/Random/RandFlat.h"

#include "QWAna/QWQVector/interface/QWQVector.h"


QWQVector::QWQVector(const edm::ParameterSet& iConfig):
	  bGen_(iConfig.getUntrackedParameter<bool>("bGen", false))
	, bSim_(iConfig.getUntrackedParameter<bool>("bSim", false))
	, bRandQ_(iConfig.getUntrackedParameter<bool>("bRandQ", false))
	, bEff_(iConfig.getUntrackedParameter<bool>("bEff", false))
	, minPt_(iConfig.getUntrackedParameter<double>("minPt", 1.0))
	, maxPt_(iConfig.getUntrackedParameter<double>("maxPt", 3.0))
	, randq_pos_(iConfig.getUntrackedParameter<double>("randq_pos", 0.5))
	, centralityToken_( iConfig.getParameter<edm::InputTag>("centrality") )
	, trackToken_(iConfig.getUntrackedParameter<edm::InputTag>("trackTag"))
	, algoParameters_(iConfig.getParameter<std::vector<int> >("algoParameters"))
	, vertexToken_( iConfig.getUntrackedParameter<edm::InputTag>("vertexSrc") )
	, epToken_( iConfig.getUntrackedParameter<edm::InputTag>("epSrc", std::string("hiEvtPlane")) )
	, fweight_( iConfig.getUntrackedParameter<edm::InputTag>("fweight", std::string("NA")) )
{

	dzdzerror_ = iConfig.getUntrackedParameter<double>("dzdzerror", 3.);
	d0d0error_ = iConfig.getUntrackedParameter<double>("d0d0error", 3.);
	pterrorpt_ = iConfig.getUntrackedParameter<double>("pterrorpt", 0.1);
	minvz_ = iConfig.getUntrackedParameter<double>("minvz", -1.);
	maxvz_ = iConfig.getUntrackedParameter<double>("maxvz", 15.);
	minEta_ = iConfig.getUntrackedParameter<double>("minEta", -2.4);
	maxEta_ = iConfig.getUntrackedParameter<double>("maxEta", 2.4);
	minCaloEta_ = iConfig.getUntrackedParameter<double>("minCaloEta", 4.4);
	maxCaloEta_ = iConfig.getUntrackedParameter<double>("maxCaloEta", 5.0);

	minCent_ = iConfig.getUntrackedParameter<int>("minCent", -1);
	maxCent_ = iConfig.getUntrackedParameter<int>("maxCent", 500);

	std::string streff = fweight_.label();
	if ( streff == std::string("NA") ) {
		std::cout << "!!! eff NA" << std::endl;
		bEff_ = false;
	} else {
		TFile * fEffFak = new TFile(streff.c_str());
		std::cout << "!!! Using particle weight " << streff << std::endl;
		if ( bEff_ ) {
			std::cout << "!!! Apply Eff correction" << std::endl;
			for ( int i = 0; i < 20; i++ ) {
				if ( streff == std::string("TrackCorrections_HIJING_538_OFFICIAL_Mar24.root") ) {
					TH2D * h = (TH2D*) fEffFak->Get("rTotalEff3D");
					for ( int c = 0; c < 200; c++ ) {
						hEff_cbin[c] = h;
					}
				}
			}
			std::cout << "!!! eff histo done" << std::endl;
		}
	}

	edm::Service<TFileService> fs;
	trV = fs->make<TTree>("trV", "trV");

	trV->Branch("pHFw", &t.pHFw, "pHFw/D");
	trV->Branch("pRe",  &t.pRe,  "pRe/D");
	trV->Branch("pIm",  &t.pIm,  "pIm/D");
	trV->Branch("nHFw", &t.nHFw, "nHFw/D");
	trV->Branch("nRe",  &t.nRe,  "nRe/D");
	trV->Branch("nIm",  &t.nIm,  "nIm/D");

	trV->Branch("pRe2", &t.pRe2,  "pRe2/D");
	trV->Branch("pIm2", &t.pIm2,  "pIm2/D");
	trV->Branch("nRe2", &t.nRe2,  "nRe2/D");
	trV->Branch("nIm2", &t.nIm2,  "nIm2/D");

	trV->Branch("cent", &(t.Cent), "cent/I");
	trV->Branch("mult", &(t.Mult), "mult/I");
	trV->Branch("noff", &(t.Noff), "mult/I");

	trV->Branch("rQ1p" , &(qval.rQ1p) , "rQ1p/D" );
	trV->Branch("rQ1p2", &(qval.rQ1p2), "rQ1p2/D");
	trV->Branch("rQ2p1", &(qval.rQ2p1), "rQ2p1/D");
	trV->Branch("rQ2p2", &(qval.rQ2p2), "rQ2p2/D");
	trV->Branch("rQ3p2", &(qval.rQ3p2), "rQ3p2/D");
	trV->Branch("rQ3p3", &(qval.rQ3p3), "rQ3p3/D");

	trV->Branch("rQ1n" , &(qval.rQ1n) , "rQ1n/D" );
	trV->Branch("rQ1n2", &(qval.rQ1n2), "rQ1n2/D");
	trV->Branch("rQ2n1", &(qval.rQ2n1), "rQ2n1/D");
	trV->Branch("rQ2n2", &(qval.rQ2n2), "rQ2n2/D");
	trV->Branch("rQ3n2", &(qval.rQ3n2), "rQ3n2/D");
	trV->Branch("rQ3n3", &(qval.rQ3n3), "rQ3n3/D");

	trV->Branch("iQ1p" , &(qval.iQ1p) , "iQ1p/D" );
	trV->Branch("iQ1p2", &(qval.iQ1p2), "iQ1p2/D");
	trV->Branch("iQ2p1", &(qval.iQ2p1), "iQ2p1/D");
	trV->Branch("iQ2p2", &(qval.iQ2p2), "iQ2p2/D");
	trV->Branch("iQ3p2", &(qval.iQ3p2), "iQ3p2/D");
	trV->Branch("iQ3p3", &(qval.iQ3p3), "iQ3p3/D");

	trV->Branch("iQ1n" , &(qval.iQ1n) , "iQ1n/D" );
	trV->Branch("iQ1n2", &(qval.iQ1n2), "iQ1n2/D");
	trV->Branch("iQ2n1", &(qval.iQ2n1), "iQ2n1/D");
	trV->Branch("iQ2n2", &(qval.iQ2n2), "iQ2n2/D");
	trV->Branch("iQ3n2", &(qval.iQ3n2), "iQ3n2/D");
	trV->Branch("iQ3n3", &(qval.iQ3n3), "iQ3n3/D");

	trV->Branch("rQMp2", &(qval.rQMp2), "rQMp2/D");
	trV->Branch("rQMn2", &(qval.rQMn2), "rQMn2/D");
	trV->Branch("iQMp2", &(qval.iQMp2), "iQMp2/D");
	trV->Branch("iQMn2", &(qval.iQMn2), "iQMn2/D");

	trV->Branch("wp1",   &(qval.wp1), "wp1/D");
	trV->Branch("wp2",   &(qval.wp2), "wp2/D");
	trV->Branch("wp3",   &(qval.wp3), "wp3/D");

	trV->Branch("wn1",   &(qval.wn1), "wn1/D");
	trV->Branch("wn2",   &(qval.wn2), "wn2/D");
	trV->Branch("wn3",   &(qval.wn3), "wn3/D");
}


//////////////////
void QWQVector::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	//std::cout << __LINE__ << "\tstart ana" << std::endl;
	if ( bGen_ ) analyzeMC(iEvent, iSetup);
	else analyzeData(iEvent, iSetup);

	if ( bSim_ ) overRide();

	std::cout << __LINE__ << "\tt.Mult = " << t.Mult << "\tt.Cent = " << t.Cent << "\tt.Noff = " << t.Noff << std::endl;
	if ( t.Mult == 0 ) return;

	edm::Handle<CaloTowerCollection> calo;
	iEvent.getByLabel("towerMaker", calo);
	t.pHFw = 0;
	t.nHFw = 0;
	t.pRe = 0;
	t.pIm = 0;
	t.pRe2 = 0;
	t.pIm2 = 0;
	t.nRe = 0;
	t.nIm = 0;
	t.nRe2 = 0;
	t.nIm2 = 0;
	if ( !calo.isValid() ) {
		std::cout << "calo failed" << std::endl;
	} else {
		for ( auto c = calo->begin(); c != calo->end(); c++ ) {
			if ( fabs(c->eta()) > maxCaloEta_ or fabs(c->eta()) < minCaloEta_ ) continue;
			if ( c->eta() > 0 ) {
				double et = c->emEt() + c->hadEt();
				t.pHFw += et;
				t.pRe += et * cos(c->phi());
				t.pIm += et * sin(c->phi());
				t.pRe2 += et * cos(2*c->phi());
				t.pIm2 += et * sin(2*c->phi());
			} else {
				double et = c->emEt() + c->hadEt();
				t.nHFw += et;
				t.nRe += et * cos(c->phi());
				t.nIm += et * sin(c->phi());
				t.nRe2 += et * cos(2*c->phi());
				t.nIm2 += et * sin(2*c->phi());
			}
		}
	}

	QHelp qh(&qval);
	qh.Fill(&t);

	trV->Fill();

	return;
}

//////////////////
QWQVector::~QWQVector()
{
	return;
}

//////////////////
void QWQVector::overRide()
{
	t.Mult = 7;
	t.Cent = 150;
	t.Charge[0] = 1;
	t.Charge[1] = 1;
	t.Charge[2] = 1;
	t.Charge[3] = -1;
	t.Charge[4] = -1;
	t.Charge[5] = -1;
	t.Charge[6] = -1;

	t.Eta[0] = 0;
	t.Eta[1] = 0;
	t.Eta[2] = 0;
	t.Eta[3] = 0;
	t.Eta[4] = 0;
	t.Eta[5] = 0;
	t.Eta[6] = 0;
	t.Eta[7] = 0;

	t.Phi[0] = 0;
	t.Phi[1] = 3.1;
	t.Phi[2] = 0.1;
	t.Phi[3] = 0.;
	t.Phi[4] = 0.1;
	t.Phi[5] = -3.1;
	t.Phi[6] = 3.1;

	t.weight[0] = 1;
	t.weight[1] = 1;
	t.weight[2] = 1;
	t.weight[3] = 1;
	t.weight[4] = 1;
	t.weight[5] = 1;
	t.weight[6] = 1;
}

///////////////////
int
QWQVector::getNoffCent(const edm::Event& iEvent, const edm::EventSetup& iSetup, int& Noff)
{
	// very hard coded Noff track centrality cut
	using namespace edm;
	using namespace reco;
//	int Noff = 0;

//	Handle<VertexCollection> vertexCollection;
//	iEvent.getByLabel(vertexToken_, vertexCollection);

	Handle<VertexCollection> vertexCollection;
	iEvent.getByLabel(vertexToken_, vertexCollection);
//	VertexCollection * recoVertices = (VertexCollection *)vertexCollection.product();
	VertexCollection recoVertices = *vertexCollection;

	if ( recoVertices.size() < 1 ) return 0;


	sort(recoVertices.begin(), recoVertices.end(), [](const reco::Vertex &a, const reco::Vertex &b){
			if ( a.tracksSize() == b.tracksSize() ) return a.chi2() < b.chi2() ? true:false;
			return a.tracksSize() > b.tracksSize() ? true:false;
			});

	int primaryvtx = 0;
	math::XYZPoint v1( (recoVertices)[primaryvtx].position().x(), (recoVertices)[primaryvtx].position().y(), (recoVertices)[primaryvtx].position().z() );
	double vxError = (recoVertices)[primaryvtx].xError();
	double vyError = (recoVertices)[primaryvtx].yError();
	double vzError = (recoVertices)[primaryvtx].zError();


	Handle<TrackCollection> tracks;
	iEvent.getByLabel(trackToken_,tracks);
	for(TrackCollection::const_iterator itTrack = tracks->begin();
		itTrack != tracks->end();
		++itTrack) {

		if ( !itTrack->quality(reco::TrackBase::highPurity) ) continue;
		if ( itTrack->charge() == 0 ) continue;
		if ( itTrack->pt() < 0.4 ) continue;

		double d0 = -1.* itTrack->dxy(v1);
		double derror=sqrt(itTrack->dxyError()*itTrack->dxyError()+vxError*vyError);
		double dz=itTrack->dz(v1);
		double dzerror=sqrt(itTrack->dzError()*itTrack->dzError()+vzError*vzError);
		if ( fabs(itTrack->eta()) > 2.4 ) continue;
		if ( fabs( dz/dzerror ) > 3. ) continue;
		if ( fabs( d0/derror ) > 3. ) continue;
		if ( itTrack->ptError()/itTrack->pt() > 0.1 ) continue;
//		bool b_pix = itTrack->numberOfValidHits() < 7;
//		if ( b_pix ) {
//			if ( fabs( dz/dzerror ) > dzdzerror_ ) continue;
//			if ( itTrack->normalizedChi2() > chi2_ ) continue;
//		} else {
//			// full track
//			if ( fabs( dz/dzerror ) > 3. ) continue;
//			if ( fabs( d0/derror ) > 3. ) continue;
//			if ( itTrack->ptError()/itTrack->pt() > pterrorpt_ ) continue;
//			if ( itTrack->numberOfValidHits() < 12 ) continue;
//		}

		Noff++;
	}

	const Int_t CentNoffCut[] = {100000, 350, 320, 300, 260, 240, 220, 185, 150, 120, 100, 80, 60, 50, 40, 30, 20, 10, 0};
	const Int_t nCentNoff = sizeof(CentNoffCut)/sizeof(Int_t);
	int cent = nCentNoff-1;
	while ( CentNoffCut[cent] <= Noff ) cent--;
	return cent;
}

//////////////////
void QWQVector::analyzeData(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	using namespace edm;
	using namespace reco;

	t.Mult = 0;
	t.Noff = 0;

	Handle<VertexCollection> vertexCollection;
	iEvent.getByLabel(vertexToken_, vertexCollection);
	VertexCollection recoVertices = *vertexCollection;
	if ( recoVertices.size() < 1 ) {
		//std::cout << __LINE__ << " vertex"<< std::endl;
		return;

	}

	sort(recoVertices.begin(), recoVertices.end(), [](const reco::Vertex &a, const reco::Vertex &b){
			if ( a.tracksSize() == b.tracksSize() ) return a.chi2() < b.chi2() ? true:false;
			return a.tracksSize() > b.tracksSize() ? true:false;
			});

	int primaryvtx = 0;
	math::XYZPoint v1( recoVertices[primaryvtx].position().x(), recoVertices[primaryvtx].position().y(), recoVertices[primaryvtx].position().z() );
	double vxError = recoVertices[primaryvtx].xError();
	double vyError = recoVertices[primaryvtx].yError();
	double vzError = recoVertices[primaryvtx].zError();

	double vz = recoVertices[primaryvtx].z();
	if (fabs(vz) < minvz_ || fabs(vz) > maxvz_) {
		//std::cout << __LINE__ << std::endl;
		return;
	}
	t.vz = vz;

	// centrality

	getNoffCent(iEvent, iSetup, t.Noff);
	if ( t.Noff < minCent_ or t.Noff >= maxCent_ ) {
		//std::cout << __LINE__ << "\tt.Noff = " << t.Noff << "\tminCent = " << minCent_ << "\tmaxCent = " << maxCent_ << std::endl;
		return;
	}
	t.Cent = t.Noff - 120;
	if ( t.Cent < 0 or t.Cent >= 200 ) {
		//std::cout << __LINE__ << std::endl;
		return;
	}

	// track
	Handle<TrackCollection> tracks;
	iEvent.getByLabel(trackToken_,tracks);

	for(TrackCollection::const_iterator itTrack = tracks->begin();
			itTrack != tracks->end();
			++itTrack) {
		if ( itTrack->charge() == 0 ) continue;
		if ( !itTrack->quality(reco::TrackBase::highPurity) ) continue;
		if ( itTrack->pt() > maxPt_ or itTrack->pt() < minPt_ ) continue;
		if ( itTrack->eta() > maxEta_ or itTrack->eta() < minEta_ ) continue;
//		if ( itTrack->numberOfValidHits() < 11 ) continue;
//		if ( itTrack->normalizedChi2() / itTrack->hitPattern().trackerLayersWithMeasurement() > 0.15 ) continue;
		if ( itTrack->ptError()/itTrack->pt() > pterrorpt_ ) continue;
		if ( itTrack->hitPattern().pixelLayersWithMeasurement() == 0 ) continue;

		double d0 = -1.* itTrack->dxy(v1);
		double derror=sqrt(itTrack->dxyError()*itTrack->dxyError()+vxError*vyError);
		if ( fabs( d0/derror ) > d0d0error_ ) continue;

		double dz=itTrack->dz(v1);
		double dzerror=sqrt(itTrack->dzError()*itTrack->dzError()+vzError*vzError);
		if ( fabs( dz/dzerror ) > dzdzerror_ ) continue;

//		if ( algoParameters_.size() != 0 and find( algoParameters_.begin(), algoParameters_.end(), itTrack->originalAlgo() ) == algoParameters_.end() ) continue;

		t.Charge[t.Mult] = itTrack->charge();
		t.Pt[t.Mult] = itTrack->pt();
		t.Eta[t.Mult] = itTrack->eta();
		t.Phi[t.Mult] = itTrack->phi();

		if ( bEff_ ) {
			double eff = hEff_cbin[t.Cent]->GetBinContent( hEff_cbin[t.Cent]->FindBin(itTrack->eta(), itTrack->pt()) );
			t.weight[t.Mult] = 1.0 / eff;
		} else {
			t.weight[t.Mult] = 1.0;
		}

		t.Mult++;
	}

	return;
}

//////////////////
void QWQVector::analyzeMC(const edm::Event&, const edm::EventSetup&)
{
	return;
}


//////////////////
void QWQVector::beginJob()
{
	return;
}

//////////////////
void QWQVector::endJob()
{
	return;
}

void QWQVector::beginRun(edm::Run const&, edm::EventSetup const&)
{
	return;
}

void QWQVector::endRun(edm::Run const&, edm::EventSetup const&)
{
	return;
}

void QWQVector::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
	return;
}

void QWQVector::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
	return;
}

//define this as a plug-in
DEFINE_FWK_MODULE(QWQVector);
