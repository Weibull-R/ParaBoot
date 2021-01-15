#include "ParaBoot.h"
#include "Pivotal.h"
#include "LSLRmodel.h"
#include <math.h>
#include <limits>

    using namespace Rcpp ;

Pivotal::Pivotal(SEXP arg1) {
// callargs<-list(ppp=x$ppp, event=event, R2=R2, CI=CI, P1=P1, P2=P2, S=S, seed=seed, dp=dp, 
// reg_order=reg_order, dist_num=dist_num, npar=npar, limit=limit,ProgRpt=ProgRpt)

	Rcpp::List L(arg1);
	positions = L["ppp"];
	event = L["event"];
	R2test = Rcpp::as<double>(L["R2"]);
	CItest = Rcpp::as<double>(L["CI"]);
	P1 = Rcpp::as<double>(L["P1"]);
	P2 = Rcpp::as<double>(L["P2"]);
	S = Rcpp::as<int>(L["S"]);
	seed = Rcpp::as<int>(L["seed"]);
	dp = L["dp"];
	regression_order = Rcpp::as<int>(L["reg_order"]);
	dist_num = Rcpp::as<int>(L["dist_num"]);
	npar = Rcpp::as<int>(L["npar"]);
	limit = Rcpp::as<double>(L["limit"]);
	ProgRpt = Rcpp::as<bool>(L["ProgRpt"]);

}

SEXP Pivotal::Execute() {
// declare the extended output objects
	int ndp = dp.size();
	arma::colvec R2(S);
	arma::mat qpiv(S,ndp);

// original PivotalMC.cpp used R2test and CItest as labels
//	because the R2 vector was so named
// variables used to control output
	int prrout=0;
	int pivout=0;
	double pvalue;
	double CCC2;
// variables in support of progrss bar, if requested
	unsigned int Spct = S/100;
	int ProgPct=0;
	int LastPct=0;
// store and return the RNG state after execution
	RNGScope scope;	
// set a consistent seed for repeatability of results
	Environment base("package:base");	
	Function SetSeed = base["set.seed"];	
	if(seed>0) {	
	SetSeed(seed);	
	}
// establish the percentiles for evaluation of confidence bounds	
	arma::colvec CBq(ndp);	
	for(int i = 0; i<ndp; i++)  {	
		CBq(i)=log(log(1/(1-dp[i])));
	}
// prepare for a recovery from taking log of values equal or less than 0
	double neginf = (-1) * std::numeric_limits<float>::max();
	arma::rowvec QCI(ndp);
	
/////////////////////////////////////////////////////////////	
/////////////// Weibull Implementation  /////////////////////
/////////////////////////////////////////////////////////////
	if(dist_num == 0)  {
		int N = event.size();
		arma::colvec y(N);
		Rcpp::NumericVector sample;
		std::vector<double> fit;
		arma::colvec coef(2);
// original pivotalMC code got away without sizing the coef vector because
// it was filled with the result of the arma::solve function.

		
		// will collect a sample covering both failures and suspensions
		// then masked by the event vector after sorting
		
		double Eta = P1;
		double Beta = P2;

		///////////// bootstrap loop  /////////
		for(unsigned int i=0; i<S; i++)  {
			// get a sample and mask it with event	
			y = Rcpp::as<arma::colvec>(rweibull(N, Beta, Eta));	
			y = arma::sort(y);	
			for(int j=0,k=0; j<N; j++)  {	
				if(event[j]==0)  {y.shed_row(j-k);}
				k++;
			}
			sample = Rcpp::wrap(y);
			
			std::unique_ptr<LSLRmodel> SAMP(new LSLRmodel(	
				sample, positions, regression_order,
				dist_num, npar, limit)
				);
			fit = SAMP->LSLRfit();	
			
			if(R2test>0.0) {			
				R2(i) = fit[npar];		// R2 position is dependent on npar
			}	

			if(CItest>0.0)  {
				if(regression_order == 0) {
				// re-defining arma coefs from Eta and Beta in the fit for 'XonY'
				coef(0) = log(fit[0]);
				coef(1) = 1/fit[1];
				}else{
				// re-defining arma coefs from Eta and Beta in the fit for 'YonX'
				coef(0) = log(fit[0])*-fit[1];
				coef(1) = fit[1];
				}

				if(npar == 2) {
				qpiv.row(i)=arma::trans((CBq-coef(0))/coef(1));
				}else{
				// how simple it was to implement 3p bounds!
				
				if(exp((CBq(0)-coef(0))/coef(1)) + fit[2] <= 0)  {
				// must build the qpiv.row(i) in a loop
				
				for(int j=0; j<ndp; j++)  {	
					if(exp((CBq(j)-coef(0))/coef(1)) + fit[2] <= 0)  {
					QCI(j) = neginf;
					}else{
					QCI(j) = log(exp((CBq(j)-coef(0))/coef(1)) + fit[2]);
					}
				}	
				qpiv.row(i)=QCI;	
				}else{	
				qpiv.row(i)=arma::trans( log(exp((CBq-coef(0))/coef(1)) + fit[2]) );	
				}
			}
			
				if(ProgRpt) {	
					// Progress report		
					ProgPct = (i+1)/Spct;		
					if(ProgPct > LastPct)  {		
						Rprintf("%3d%% completion",(i+1)/Spct);		
						Rprintf("\r");		
						R_FlushConsole();		
						R_ProcessEvents();		
					}		
					LastPct = ProgPct;		
				}
			}
		//close bootstrap loop		
		}		
			
	// close dist_num == 0		
	}		

	int LCB=0;
	arma::rowvec LBpiv(ndp);
	arma::rowvec HBpiv(ndp);
	arma::rowvec median(ndp);
	
// process the prr vector according to output control					
	if(R2test>0.0) {				
	prrout = 1;				
	R2=arma::sort(R2);				
	if(R2test< 1.0) {				
		prrout=2;			
		arma::colvec Absolute(1);			
		Absolute(0)=1.0;			
		R2=join_cols(R2,Absolute);			
// pve_u is an integer representation  of the percentile of R-square					
		arma::uvec pvalue_u=arma::find(R2>R2test,1,"first");			
// as long as S is sufficiently large (>10^4) there is no accuracy to be gained by interpolation					
		pvalue= (double) (pvalue_u(0)) /S*100;			
// note: integer math in this dimension specification may breakdown if S!= Multiple of 10					
		CCC2= (double) R2(S/10-1);			
					
					
	}}				
					
	if(CItest>0.0) {				
		pivout=1;			
		for(int i=0; i<ndp; i++)  {			
			qpiv.col(i)=arma::sort(qpiv.col(i));		
		}			
		if(CItest< 1.0) {				
			pivout=2;			
			LCB=(int) S*(1-CItest)/2;			
			LBpiv=qpiv.row(LCB-1);			
			HBpiv=qpiv.row(S-LCB-1);			
			median=qpiv.row(S/2-1);
		}					
	}				
					
	int outputcase=prrout+4*pivout;				
	switch(outputcase)				
	{				
		case 1:			
			return wrap(R2);		
			break;
			
		case 2:			
			return DataFrame::create(			
				Rcpp::Named("Pvalue")=pvalue,		
				Rcpp::Named("CCC2")=CCC2 );		
			break;		
					
		case 4:				
			return wrap(qpiv);			
			break;			
	// this is the unlikely case that both extended output objects are called for					
		case 5:				
				return List::create(			
				Rcpp::Named("prr")=wrap(R2),			
				Rcpp::Named("pivotals")=wrap(qpiv) );			
			break;			
						
		case 6:				
			return List::create(			
				Rcpp::Named("prrCCC2")=			
				DataFrame::create(			
					Rcpp::Named("Pvalue")=pvalue,		
					Rcpp::Named("CCC2")=CCC2 ),		
				Rcpp::Named("pivotals")=wrap(qpiv) );			
			break;			
						
		case 8:				
			return DataFrame::create(			
				Rcpp::Named("Lower")=wrap(arma::trans(LBpiv)),		
				Rcpp::Named("Median")=wrap(arma::trans(median)),		
				Rcpp::Named("Upper")=wrap(arma::trans(HBpiv)) );		
			break;			
						
		case 9:				
			return List::create(			
				Rcpp::Named("prr")=wrap(R2),		
				Rcpp::Named("pivotals")=		
				DataFrame::create(		
					Rcpp::Named("Lower")=wrap(arma::trans(LBpiv)),	
					Rcpp::Named("Median")=wrap(arma::trans(median)),	
					Rcpp::Named("Upper")=wrap(arma::trans(HBpiv)) )	
				);		
			break;			
						
		case 10:				
			return List::create(			
				Rcpp::Named("prrCCC2")=		
				 DataFrame::create(		
					Rcpp::Named("Pvalue")=pvalue,	
					Rcpp::Named("CCC2")=CCC2 ),	
				Rcpp::Named("pivotals")=		
				DataFrame::create(		
					Rcpp::Named("Lower")=wrap(arma::trans(LBpiv)),	
					Rcpp::Named("Median")=wrap(arma::trans(median)),	
					Rcpp::Named("Upper")=wrap(arma::trans(HBpiv)) )	
				);		
			break;			
						
		default:			
			return wrap(0.0);			
					
	}				
				
}

// external functions
SEXP pivotalMC2(SEXP arg1) {
// Construct an LSLRmodel with this input argument
	std::unique_ptr<Pivotal> PIV(new Pivotal(arg1));
	return PIV->Execute();
}

