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
	, bPPreco_(iConfig.getUntrackedParameter<bool>("bPPreco", false))
	, bRandQ_(iConfig.getUntrackedParameter<bool>("bRandQ", false))
	, bEff_(iConfig.getUntrackedParameter<bool>("bEff", false))
	, minPt_(iConfig.getUntrackedParameter<double>("minPt", 1.0))
	, maxPt_(iConfig.getUntrackedParameter<double>("maxPt", 3.0))
	, randq_pos_(iConfig.getUntrackedParameter<double>("randq_pos", 0.5))
	, centralityToken_( consumes<int>(iConfig.getParameter<edm::InputTag>("centrality")) )
	, trackToken_(consumes<reco::TrackCollection>(iConfig.getUntrackedParameter<edm::InputTag>("trackTag")))
	, algoParameters_(iConfig.getParameter<std::vector<int> >("algoParameters"))
	, vertexToken_( consumes<reco::VertexCollection>(iConfig.getUntrackedParameter<edm::InputTag>("vertexSrc")) )
	, epToken_( consumes<reco::EvtPlaneCollection>(iConfig.getUntrackedParameter<edm::InputTag>("epSrc", std::string("hiEvtPlane") )) )
	, fweight_( iConfig.getUntrackedParameter<edm::InputTag>("fweight", std::string("NA")) )
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

	if ( bPPreco_ ) {
		towerToken_ = consumes<CaloTowerCollection>(iConfig.getUntrackedParameter<edm::InputTag>("towerMaker"));
	}

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
				if ( streff == std::string("Hydjet_eff_mult_v1.root") ) {
					TH2D * h = (TH2D*) fEffFak->Get("rTotalEff3D_1");
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

	//std::cout << __LINE__ << "\t" << t.Mult << std::endl;
	if ( t.Mult == 0 ) return;

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

//////////////////
//int QWQVector::getNoffCent(const edm::Event& iEvent, const edm::EventSetup& iSetup, int& Noff)
//{
//	// very hard coded Noff track centrality cut
//	using namespace edm;
//	using namespace reco;
//	//      int Noff = 0;
//
//	Handle<VertexCollection> vertexCollection;
//	iEvent.getByToken(vertexToken_, vertexCollection);
//	const VertexCollection * recoVertices = vertexCollection.product();
//
//	if ( recoVertices.size() < 1 ) return;
//	sort(recoVertices.begin(), recoVertices.end(), [](const reco::Vertex &a, const reco::Vertex &b){
//			if ( a.tracksSize() == b.tracksSize() ) return a.chi2() < b.chi2() ? true:false;
//			return a.tracksSize() > b.tracksSize() ? true:false;
//			});
//
//	int primaryvtx = 0;
//	math::XYZPoint v1( (*recoVertices)[primaryvtx].position().x(), (*recoVertices)[primaryvtx].position().y(), (*recoVertices)[primaryvtx].position().z() );
//	double vxError = (*recoVertices)[primaryvtx].xError();
//	double vyError = (*recoVertices)[primaryvtx].yError();
//	double vzError = (*recoVertices)[primaryvtx].zError();
//
//
//	Handle<TrackCollection> tracks;
//	iEvent.getByToken(trackToken_,tracks);
//	for(TrackCollection::const_iterator itTrack = tracks->begin();
//			itTrack != tracks->end();
//			++itTrack) {
//
//		if ( !itTrack->quality(reco::TrackBase::highPurity) ) continue;
//		if ( itTrack->charge() == 0 ) continue;
//		if ( itTrack->pt() < 0.4 ) continue;
//
//		double d0 = -1.* itTrack->dxy(v1);
//		double derror=sqrt(itTrack->dxyError()*itTrack->dxyError()+vxError*vyError);
//		double dz=itTrack->dz(v1);
//		double dzerror=sqrt(itTrack->dzError()*itTrack->dzError()+vzError*vzError);
//		if ( fabs(itTrack->eta()) > 2.4 ) continue;
//		if ( fabs( dz/dzerror ) > 3. ) continue;
//		if ( fabs( d0/derror ) > 3. ) continue;
//		if ( itTrack->ptError()/itTrack->pt() > 0.1 ) continue;
//
//		Noff++;
//	}
//
//	int cent = ;
//	while ( [cent] <= Noff ) cent--;
//	return cent;
//}
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
		//std::cout << __LINE__ << std::endl;
		return;
	}
	t.vz = vz;

	// centrality

	if ( bPPreco_ ) {
		double etHFtowerSumPlus=0;
		double etHFtowerSumMinus=0;
		double etHFtowerSum=0;
		Handle<CaloTowerCollection> towers;
		iEvent.getByToken(towerToken_,towers);
		for( size_t i = 0; i<towers->size(); ++ i){
			const CaloTower & tower = (*towers)[ i ];
			double eta = tower.eta();
			bool isHF = tower.ietaAbs() > 29;
			if(isHF && eta > 0){
				etHFtowerSumPlus += tower.pt();
			}
			if(isHF && eta < 0){
				etHFtowerSumMinus += tower.pt();
			}
		}
		etHFtowerSum=etHFtowerSumPlus + etHFtowerSumMinus;

		double binLowEdge[200]={4487.37, 4370.52, 4279.22, 4199.03, 4113.36, 4030.29, 3947.37, 3870.82, 3797.94, 3721.48, 3648.81, 3575.99, 3506.18, 3440.74, 3374.5, 3310.49, 3249.72, 3190.49, 3127.43, 3066.91, 3012.16, 2954.08, 2897.16, 2840.3, 2786.54, 2735.06, 2682.83, 2631.95, 2580.71, 2529.93, 2483.34, 2436.59, 2389.05,  2343.58, 2300.27, 2256.49, 2210.35, 2167.14, 2128.09, 2086.24, 2044.85, 2002.72, 1962.42, 1925.23, 1889.2, 1851.68, 1815.58, 1778.47, 1743.48, 1706.47, 1671.08,  1636.7, 1604.94, 1571.63, 1539.86, 1508.37, 1477.12, 1445.73, 1417.7, 1387.98, 1359.02, 1330.3, 1301.45, 1274.07, 1246.54, 1219.36, 1191.97, 1165.77, 1140.4,     1114.92, 1091.98, 1067.94, 1043.67, 1019.66, 995.39, 970.466, 947.786, 924.75, 902.723, 879.824, 859.262, 838.212, 817.18, 796.627, 776.494, 757.142, 737.504,    719.604, 701.142, 684.043, 665.89, 648.427, 630.224, 612.877, 596.435, 580.397, 565.396, 550.272, 535.204, 520.48, 505.854, 491.648, 477.531, 463.192, 449.773,   436.806, 423.944, 410.4, 397.962, 386.135, 374.47, 362.499, 351.17, 339.635, 328.402, 317.875, 307.348, 296.957, 287.002, 276.94, 267.822, 258.796, 249.366, 239.974, 231.563, 223.362, 214.902, 206.818, 199.417, 191.609, 184.184, 177.042, 169.839, 163.579, 157.186, 151.136, 145.165, 139.213, 133.218, 127.748, 122.445, 117.458, 112.715, 108.179, 103.713, 99.2518, 94.8864, 90.7892, 86.692, 82.819, 79.0331, 75.4791, 71.8774, 68.5738, 65.5363, 62.6369, 59.7441, 57.0627, 54.3838, 51.7242, 49.1577, 46.7914, 44.4615, 42.3374, 40.2863, 38.2674, 36.3979, 34.4769, 32.7274, 30.9911, 29.3998, 27.7739, 26.2442, 24.795, 23.3496, 21.8717, 20.5263, 19.2405, 18.08, 16.9542, 15.882, 14.8344, 13.8014, 12.7824, 11.8165, 10.8308, 9.94351, 9.08363, 8.20773, 7.40535, 6.57059, 5.81859, 5.0626, 4.32634, 3.57026, 2.83467, 2.09189, 1.36834, 0.673038, 0};

		t.Cent = -1;
		for(int i=0; i<200; i++){
			if(etHFtowerSum>=binLowEdge[i]){
				t.Cent = i; break;
			}
		}
	} else {
		edm::Handle<int> ch;
		iEvent.getByToken(centralityToken_,ch);
		t.Cent = *(ch.product());
		if ( t.Cent < 0 or t.Cent >= 200 ) {
			return;
		}
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

		if ( algoParameters_.size() != 0 and find( algoParameters_.begin(), algoParameters_.end(), itTrack->originalAlgo() ) == algoParameters_.end() ) continue;

		if ( bRandQ_ ) {
			edm::Service<edm::RandomNumberGenerator> rng;
			CLHEP::HepRandomEngine* engine = &rng->getEngine(iEvent.streamID());
			if ( CLHEP::RandFlat::shoot(engine) < randq_pos_ ) {
				t.Charge[t.Mult] = 1;
			} else {
				t.Charge[t.Mult] = -1;
			}
		} else {
			t.Charge[t.Mult] = itTrack->charge();
		}
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
