
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


#include "QWAna/QWQVector/interface/QWQVector.h"


QWQVector::QWQVector(const edm::ParameterSet& iConfig):
	  bGen_(iConfig.getUntrackedParameter<bool>("bGen", false))
	, minPt_(iConfig.getUntrackedParameter<bool>("minPt", 1.0))
	, maxPt_(iConfig.getUntrackedParameter<bool>("maxPt", 3.0))
	, trackTag_(iConfig.getUntrackedParameter<edm::InputTag>("trackTag"))
	, algoParameters_(iConfig.getParameter<std::vector<int> >("algoParameters"))
	, vertexToken_( consumes<reco::VertexCollection>(iConfig.getUntrackedParameter<edm::InputTag>("vertexSrc")) )
	, epToken_( consumes<reco::EvtPlaneCollection>(iConfig.getUntrackedParameter<edm::InputTag>("epSrc", string("hiEvtPlane") )) )
{
	trackToken_ = consumes<reco::TrackCollection>(trackTag_);

	dzdzerror_ = iConfig.getUntrackedParameter<double>("dzdzerror", 3.);
	d0d0error_ = iConfig.getUntrackedParameter<double>("d0d0error", 3.);
	pterrorpt_ = iConfig.getUntrackedParameter<double>("pterrorpt", 0.1);
	minvz_ = iConfig.getUntrackedParameter<double>("minvz", -1.);
	maxvz_ = iConfig.getUntrackedParameter<double>("maxvz", 15.);

	minCent_ = iConfig.getUntrackedParameter<int>("minCent", -1);
	maxCent_ = iConfig.getUntrackedParameter<int>("maxCent", 500);


	edm::Service<TFileService> fs;
	trV = fs->make<TTree>("trV", "trV");

	trV->Branch("pw",  pw,  "pw[12]/D");
	trV->Branch("pre", pre, "pre[12]/D");
	trV->Branch("pim", pim, "pim[12]/D");

	trV->Branch("nw",  nw,  "nw[12]/D");
	trV->Branch("nre", nre, "nre[12]/D");
	trV->Branch("nim", nim, "nim[12]/D");

	trV->Branch("pHFw", &t.pHFw, "pHFw/D");
	trV->Branch("pRe",  &t.pRe,  "pRe/D");
	trV->Branch("pIm",  &t.pIm,  "pIm/D");

	trV->Branch("nHFw", &t.nHFw, "nHFw/D");
	trV->Branch("nRe",  &t.nRe,  "nRe/D");
	trV->Branch("nIm",  &t.nIm,  "nIm/D");

	trV->Branch("cent", &(t.Cent), "cent/I");

}


//////////////////
void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	if ( bGen_ ) analyzeMC(iEvent, iSetup);
	else analyzeData(iEvent, iSetup);

	if ( t.Mult == 0 ) return;

	QVector pq[12];
	QVector nq[12];

	for ( int i = 0; i < t.Mult; i++ ) {
		int bin = (t.Eta[i] + 2.4) / 0.4;
		if ( bin < 0 or bin > 11 ) continue;
		if ( t.Charge[i] > 0 ) {
			pq[bin].AddParticle(t.Phi[i], t.Pt[i]);
		} else {
			nq[bin].AddParticle(t.Phi[i], t.Pt[i]);
		}
	}

	for ( int i = 0; i < 12; i++ ) {
		pw[i] =  pq[i].GetW();
		pre[i] = pq[i].GetQ().real();
		pim[i] = pq[i].GetQ().imag();
		nw[i] =  nq[i].GetW();
		nre[i] = nq[i].GetQ().real();
		nim[i] = nq[i].GetQ().imag();
	}
	trV->Fill();

	return;
}


//////////////////
void analyzeData(const edm::Event&, const edm::EventSetup&)
{
	using namespace edm;
	using namespace reco;

	t.Mult = 0;

	Handle<VertexCollection> vertexCollection;
	iEvent.getByToken(vertexToken_, vertexCollection);
	VertexCollection recoVertices = *vertexCollection;
	if ( recoVertices.size() > nvtx_ ) return;
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
		return;
	}
	t.vz = vz;

	// centrality
	edm::Handle<int> ch;
	iEvent.getByToken(centralityToken_,ch);
	t.Cent = *(ch.product());
	if ( t.Cent < 0 or t.Cent >= 200 ) {
		return;
	}

	// track
	Handle<TrackCollection> tracks;
	iEvent.getByToken(trackToken_,tracks);

	for(TrackCollection::const_iterator itTrack = tracks->begin();
			itTrack != tracks->end();
			++itTrack) {
		if ( itTrack->charge() == 0 ) continue;
		if ( !itTrack->quality(reco::TrackBase::highPurity) ) continue;
		if ( itTrack->pt() > maxPt_ or itTrack->pt() < minPt_ ) continue;
		if ( fabs(itTrack->eta()) > 2.4 ) continue;
		if ( itTrack->numberOfValidHits() < 11 ) continue;
		if ( itTrack->normalizedChi2() / itTrack->hitPattern().trackerLayersWithMeasurement() > 0.15 ) continue;
		if ( itTrack->ptError()/itTrack->pt() > pterrorpt_ ) continue;

		double d0 = -1.* itTrack->dxy(v1);
		double derror=sqrt(itTrack->dxyError()*itTrack->dxyError()+vxError*vyError);
		if ( fabs( d0/derror ) > d0d0error_ ) continue;

		double dz=itTrack->dz(v1);
		double dzerror=sqrt(itTrack->dzError()*itTrack->dzError()+vzError*vzError);
		if ( fabs( dz/dzerror ) > dzdzerror_ ) continue;

		if ( find( algoParameters_.begin(), algoParameters_.end(), itTrack->originalAlgo() ) == algoParameters_.end() ) continue;

		t->Charge[t->Mult] = itTrack->charge();
		t->Pt[t->Mult] = itTrack->pt();
		t->Eta[t->Mult] = itTrack->eta();
		t->Phi[t->Mult] = itTrack->phi();
		t->Mult++;
	}

        edm::Handle<reco::EvtPlaneCollection> epCollection;
        iEvent.getByToken(epToken_, epCollection);
        if ( ! epCollection.isValid() ) return;
        const reco::EvtPlaneCollection * ep = epCollection.product();
        if ( ep->size() != hi::NumEPNames ) return;
	t.nRe = (*ep)[hi::HFm2].qx(0);
	t.nIm = (*ep)[hi::HFm2].qy(0);
	t.nHFw = (*ep)[hi::HFm2].sumw();

	t.pRe = (*ep)[hi::HFp2].qx(0);
	t.pIm = (*ep)[hi::HFp2].qy(0);
	t.pHFw = (*ep)[hi::HFp2].sumw();

	return;
}

//////////////////
void analyzeMC(const edm::Event&, const edm::EventSetup&)
{

	return;
}

//////////////////
void beginJob()
{
	return;
}

//////////////////
void endJob()
{
	return;
}

//////////////////
void beginRun(edm::Run const&, edm::EventSetup const&)
{
	return;
}

//////////////////
void endRun(edm::Run const&, edm::EventSetup const&)
{
	return;
}

//////////////////
void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
	return;
}

//////////////////
void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
	return;
}


