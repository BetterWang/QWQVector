#include <complex>
#include "TTree.h"

// Q Vector helper class
class QVector {
	typedef std::complex<double> Complex;
public:
	QVector():q_(0,0), weight_(0){};
	QVector(double x, double y, double w):q_(x,y), weight_(w){};
	void AddParticle(double phi, double w) {
		q_ += w * Complex(cos(phi), sin(phi));
		weight_ += w;
	};
	void RemoveParticle(double phi, double w) {
		q_ -= w * Complex(cos(phi), sin(phi));
		weight_ -= w;
	};
	Complex GetQ() {return q_;};
	double GetW() {return weight_;};

private:
	Complex q_;
	double weight_;
};

// event structure
const int NMAX_TRK = 10000;
typedef struct QWEvent_ {
        int     Cent;
        int     Mult;
        double  vz;
        int     Noff;
        double  Pt[NMAX_TRK];
        double  Eta[NMAX_TRK];
        double  Phi[NMAX_TRK];
        int     Charge[NMAX_TRK];
        double  weight[NMAX_TRK];
        int     RunId;
        int     EventId;
	double	pRe;
	double	pIm;
	double	nRe;
	double	nIm;
	double	pHFw;
	double	nHFw;
} QWEvent;

class QWQVector : public edm::EDAnalyzer {
public:
	explicit QWQVector(const edm::ParameterSet&);
	~QWQVector();

private:
	virtual void beginJob();
	virtual void analyze(const edm::Event&, const edm::EventSetup&);
	virtual void endJob();

	virtual void beginRun(edm::Run const&, edm::EventSetup const&);
	virtual void endRun(edm::Run const&, edm::EventSetup const&);
	virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
	virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

	////////////////////////////////
	void analyzeData(const edm::Event&, const edm::EventSetup&);
	void analyzeMC(const edm::Event&, const edm::EventSetup&);


	bool bGen_;

	double minPt_;
	double maxPt_;

	edm::EDGetTokenT<int>                           centralityToken_;
	edm::EDGetTokenT<reco::TrackCollection>		trackToken_;
	std::vector<int>				algoParameters_;
	edm::EDGetTokenT<reco::VertexCollection>	vertexToken_;
	edm::EDGetTokenT<reco::EvtPlaneCollection>      epToken_;


	double	dzdzerror_;
	double	d0d0error_;
	double	pterrorpt_;
	double  minvz_;
	double  maxvz_;
	int	minCent_;
	int	maxCent_;



	QWEvent t;

	TTree * trV;

	double pw[12];
	double pre[12];
	double pim[12];
	double nw[12];
	double nre[12];
	double nim[12];

};
