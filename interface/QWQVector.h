#include <correlations/Types.hh>
#include <correlations/Result.hh>
#include <correlations/QVector.hh>
#include <correlations/closed/FromQVector.hh>

#include <complex>
#include <utility>
#include "TTree.h"

typedef std::complex<double> Complex;
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

	double	pRe2;
	double	pIm2;
	double	nRe2;
	double	nIm2;
} QWEvent;

// Q Vector helper class

typedef struct QValue_ {

	double iQ1p,  iQ1n;
	double iQ2p1, iQ2p2;
	double iQ2n1, iQ2n2;
	double iQ3p2, iQ3p3;
	double iQ3n2, iQ3n3;

	double rQ1p,  rQ1n;
	double rQ2p1, rQ2p2;
	double rQ2n1, rQ2n2;
	double rQ3p2, rQ3p3;
	double rQ3n2, rQ3n3;

	double wp1, wn1;
	double wp2, wn2;
	double wp3, wn3;
} QValue;

class QHelp {
public:
	QHelp( QValue* q ) : qval(q), h1{1}, h2{2, -2}, h3{1, 1, -2},
	       q1p1(h1, true), q2p2(h2, true), q3p3(h3, true),
	       q1n1(h1, true), q2n2(h2, true), q3n3(h3, true) {};
	void Fill(QWEvent * const t) {
		for ( int i = 0; i < t->Mult; i++ ) {
			if ( t->Charge[i] > 0 ) {
				q1p1.fill(t->Phi[i], 1.);
				q2p2.fill(t->Phi[i], 1.);
				q3p3.fill(t->Phi[i], 1.);
			} else {
				q1n1.fill(t->Phi[i], 1.);
				q2n2.fill(t->Phi[i], 1.);
				q3n3.fill(t->Phi[i], 1.);
			}
		}
		GetQ();
	};

	void GetQ() {
		correlations::closed::FromQVector cq1p1(q1p1);
		correlations::closed::FromQVector cq2p2(q2p2);
		correlations::closed::FromQVector cq3p3(q3p3);
		correlations::closed::FromQVector cq1n1(q1n1);
		correlations::closed::FromQVector cq2n2(q2n2);
		correlations::closed::FromQVector cq3n3(q3n3);

		correlations::Result r1p1 = cq1p1.calculate(1, h1);
		correlations::Result r1n1 = cq1n1.calculate(1, h1);

		correlations::Result r2p1 = cq2p2.calculate(1, h2);
		correlations::Result r2p2 = cq2p2.calculate(2, h2);
		correlations::Result r2n1 = cq2n2.calculate(1, h2);
		correlations::Result r2n2 = cq2n2.calculate(2, h2);

		correlations::Result r3p2 = cq3p3.calculate(2, h3);
		correlations::Result r3p3 = cq3p3.calculate(3, h3);
		correlations::Result r3n2 = cq3n3.calculate(2, h3);
		correlations::Result r3n3 = cq3n3.calculate(3, h3);

		qval->rQ1p  = r1p1.sum().real();
		qval->rQ1n  = r1n1.sum().real();

		qval->rQ2p1 = r2p1.sum().real();
		qval->rQ2p2 = r2p2.sum().real();
		qval->rQ2n1 = r2n1.sum().real();
		qval->rQ2n2 = r2n2.sum().real();

		qval->rQ3p2 = r3p2.sum().real();
		qval->rQ3p3 = r3p3.sum().real();
		qval->rQ3n2 = r3n2.sum().real();
		qval->rQ3n3 = r3n3.sum().real();

		qval->iQ1p  = r1p1.sum().imag();
		qval->iQ1n  = r1n1.sum().imag();

		qval->iQ2p1 = r2p1.sum().imag();
		qval->iQ2p2 = r2p2.sum().imag();
		qval->iQ2n1 = r2n1.sum().imag();
		qval->iQ2n2 = r2n2.sum().imag();

		qval->iQ3p2 = r3p2.sum().imag();
		qval->iQ3p3 = r3p3.sum().imag();
		qval->iQ3n2 = r3n2.sum().imag();
		qval->iQ3n3 = r3n3.sum().imag();


		qval->wp1 = r1p1.weight();
		qval->wn1 = r1n1.weight();
		qval->wp2 = r2p2.weight();
		qval->wn2 = r2n2.weight();
		qval->wp3 = r3p3.weight();
		qval->wn3 = r3n3.weight();

	};

private:
	QValue* qval;
	correlations::HarmonicVector h1;
	correlations::HarmonicVector h2;
	correlations::HarmonicVector h3;
	correlations::QVector	q1p1, q2p2, q3p3;
	correlations::QVector	q1n1, q2n2, q3n3;
};

class QVector {
public:
	QVector(int N):N_(N), q_(0,0), weight_(0){};
	QVector(int N, double x, double y, double w):N_(N), q_(x,y), weight_(w){};
	void AddParticle(double phi, double w=1.) {
		q_ += w * Complex(cos(N_*phi), sin(N_*phi));
		weight_ += w;
	};
	void RemoveParticle(double phi, double w) {
		q_ -= w * Complex(cos(N_*phi), sin(N_*phi));
		weight_ -= w;
	};
	Complex GetQ() {return q_;};
	double GetW() {return weight_;};
	int GetN() { return N_; };

private:
	int N_;
	Complex q_;
	double weight_;
};

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
	double  minEta_;
	double  maxEta_;

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

	double pre2[12];
	double pim2[12];
	double nre2[12];
	double nim2[12];

	QValue qval;

};
