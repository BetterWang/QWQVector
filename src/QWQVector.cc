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


#include "QWAna/QWQVector/interface/QWQVector.h"


QWQVector::QWQVector(const edm::ParameterSet& iConfig):
	  bGen_(iConfig.getUntrackedParameter<bool>("bGen", false))
	, minPt_(iConfig.getUntrackedParameter<double>("minPt", 1.0))
	, maxPt_(iConfig.getUntrackedParameter<double>("maxPt", 3.0))
	, centralityToken_( consumes<int>(iConfig.getParameter<edm::InputTag>("centrality")) )
	, trackToken_(consumes<reco::TrackCollection>(iConfig.getUntrackedParameter<edm::InputTag>("trackTag")))
	, algoParameters_(iConfig.getParameter<std::vector<int> >("algoParameters"))
	, vertexToken_( consumes<reco::VertexCollection>(iConfig.getUntrackedParameter<edm::InputTag>("vertexSrc")) )
	, epToken_( consumes<reco::EvtPlaneCollection>(iConfig.getUntrackedParameter<edm::InputTag>("epSrc", std::string("hiEvtPlane") )) )
{

	dzdzerror_ = iConfig.getUntrackedParameter<double>("dzdzerror", 3.);
	d0d0error_ = iConfig.getUntrackedParameter<double>("d0d0error", 3.);
	pterrorpt_ = iConfig.getUntrackedParameter<double>("pterrorpt", 0.1);
	minvz_ = iConfig.getUntrackedParameter<double>("minvz", -1.);
	maxvz_ = iConfig.getUntrackedParameter<double>("maxvz", 15.);
	minEta_ = iConfig.getUntrackedParameter<double>("minEta", -2.4);
	maxEta_ = iConfig.getUntrackedParameter<double>("maxEta", 2.4);

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

	trV->Branch("pre2", pre, "pre2[12]/D");
	trV->Branch("pim2", pim, "pim2[12]/D");
	trV->Branch("nre2", nre, "nre2[12]/D");
	trV->Branch("nim2", nim, "nim2[12]/D");

	trV->Branch("pHFw", &t.pHFw, "pHFw/D");
	trV->Branch("pRe",  &t.pRe,  "pRe/D");
	trV->Branch("pIm",  &t.pIm,  "pIm/D");
	trV->Branch("nHFw", &t.nHFw, "nHFw/D");
	trV->Branch("nRe",  &t.nRe,  "nRe/D");
	trV->Branch("nIm",  &t.nIm,  "nIm/D");

	trV->Branch("pRe2",  &t.pRe2,  "pRe2/D");
	trV->Branch("pIm2",  &t.pIm2,  "pIm2/D");
	trV->Branch("nRe2",  &t.nRe2,  "nRe2/D");
	trV->Branch("nIm2",  &t.nIm2,  "nIm2/D");

	trV->Branch("cent", &(t.Cent), "cent/I");
	trV->Branch("mult", &(t.Mult), "mult/I");

	qval = new QValue();
	trV->Branch("QValue", qval, 8000, 1);
}


//////////////////
void QWQVector::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	if ( bGen_ ) analyzeMC(iEvent, iSetup);
	else analyzeData(iEvent, iSetup);


	if ( t.Mult == 0 ) return;

	std::vector<QVector> pq(12, QVector(1));
	std::vector<QVector> nq(12, QVector(1));
	std::vector<QVector> pq2(12, QVector(2));
	std::vector<QVector> nq2(12, QVector(2));

	for ( int i = 0; i < t.Mult; i++ ) {
		int bin = (t.Eta[i] + 2.4) / 0.4;
		if ( bin < 0 or bin > 11 ) continue;
		if ( t.Charge[i] > 0 ) {
			pq[bin].AddParticle(t.Phi[i]);
			pq2[bin].AddParticle(t.Phi[i]);
		} else {
			nq[bin].AddParticle(t.Phi[i]);
			nq2[bin].AddParticle(t.Phi[i]);
		}
	}

	for ( int i = 0; i < 12; i++ ) {
		pw[i] =  pq[i].GetW();
		pre[i] = pq[i].GetQ().real();
		pim[i] = pq[i].GetQ().imag();
		nw[i] =  nq[i].GetW();
		nre[i] = nq[i].GetQ().real();
		nim[i] = nq[i].GetQ().imag();

		pre2[i] = pq2[i].GetQ().real();
		pim2[i] = pq2[i].GetQ().imag();
		nre2[i] = nq2[i].GetQ().real();
		nim2[i] = nq2[i].GetQ().imag();
	}

	QHelp qh(qval);

	qh.Fill( &t );
	qh.GetQ();

	trV->Fill();

	return;
}

//////////////////
QWQVector::~QWQVector()
{
	return;
}

//////////////////
void QWQVector::analyzeData(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	using namespace edm;
	using namespace reco;

	t.Mult = 0;

	Handle<VertexCollection> vertexCollection;
	iEvent.getByToken(vertexToken_, vertexCollection);
	VertexCollection recoVertices = *vertexCollection;
	if ( recoVertices.size() < 1 ) return;
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

		t.Charge[t.Mult] = itTrack->charge();
		t.Pt[t.Mult] = itTrack->pt();
		t.Eta[t.Mult] = itTrack->eta();
		t.Phi[t.Mult] = itTrack->phi();
		t.Mult++;
	}

	edm::Handle<reco::EvtPlaneCollection> epCollection;
	iEvent.getByToken(epToken_, epCollection);
	if ( ! epCollection.isValid() ) return;
	const reco::EvtPlaneCollection * ep = epCollection.product();
	if ( ep->size() != hi::NumEPNames ) return;

	t.nRe2 = (*ep)[hi::HFm2].qx(0);
	t.nIm2 = (*ep)[hi::HFm2].qy(0);
	t.pRe2 = (*ep)[hi::HFp2].qx(0);
	t.pIm2 = (*ep)[hi::HFp2].qy(0);

	t.nRe = (*ep)[hi::HFm1].qx(0);
	t.nIm = (*ep)[hi::HFm1].qy(0);
	t.nHFw = (*ep)[hi::HFm1].sumw();
	t.pRe = (*ep)[hi::HFp1].qx(0);
	t.pIm = (*ep)[hi::HFp1].qy(0);
	t.pHFw = (*ep)[hi::HFp1].sumw();

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

