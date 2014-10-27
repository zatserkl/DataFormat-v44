// Class Wepl members are constructor Wepl,and function EtoWEPL 
// converting  energy deposited in stages (in MeV) to WEPL in mm.
// Returned WEPL equal to -1000   means that energy deposition
// pattern is inconsistent with Bragg curve for stages 1-4. 
// Returned WEPL equal to 1000 means that energy deposition in a  
// stage is well above maximum for a single proton event.
// WEPL = 2000 if enrgy depositions in all stages are below thresholds.
// Function SetEthresholds is used to set energy thresholds for 5 stages
// and auxiliary function myP4 calculates 4 degree polinomial.
// Created by V.A.Bashkirov, May 9,2014, mdified May 16, 2014
// Jul 24, 2014. Tuned for stage interfaces handling  for Jul 20 beamtest data
 class Wepl {  
   Float_t p0[5],p0c[5],p1[5],p1c[5],p2[5],p2c[5],p3[5],p3c[5],p4[5];
	Float_t thr0,thr1,thr2,thr3,thr4;	  

  public: 
	Wepl(Float_t * );		//constructor
	Float_t EtoWEPL(Float_t *);  // returns WEPL in mm
	Float_t myP4(Float_t *, Float_t);
	void SetEthresholds(Float_t,Float_t,Float_t,Float_t,Float_t);
	
};
Wepl::Wepl(Float_t * Fpar ) {
	for (int i=0; i<5; ++i)  {
	p0[i]=Fpar[i]; p0c[i]=Fpar[i+5];p1[i]=Fpar[i+10];p1c[i]=Fpar[i+15];
	p2[i]=Fpar[i+20]; p2c[i]=Fpar[i+25];p3[i]=Fpar[i+30];p3c[i]=Fpar[i+35];
	p4[i]=Fpar[i+40];        }
}
Float_t Wepl::EtoWEPL(Float_t* Estage) {
	Float_t wet_wm,wet_pr;	
	if(Estage[4]>thr4) { if(Estage[4]>87.) return 1000;
 	wet_wm=myP4(p4,Estage[4]); 
 //Check if energy deposition pattern is compatible to Bragg curve (within 5 sigma): 	
	if(Estage[3]<35 || Estage[2]<25 || Estage[1]<21 || Estage[0]<19)  return -1000;
 	return wet_wm*1.0383; // polystyrene to water equivalent
		}
	if(Estage[3]>thr3) { if(Estage[3]>87.) return 1000;
 	wet_wm=myP4(p3,Estage[3]);
	if(Estage[2]<25 || Estage[1]<21 || Estage[0]<19)  return -1000;
 	if(Estage[3]>69.5) { //If Bragg peak shared by stage3,wrapping, and Estage[4]<trh4
 	wet_pr=myP4(p2c,Estage[2]); 
 	wet_wm=wet_wm*.3+wet_pr*.7; //use weighted mean with sigma 4 and 6mm
 		     }   
	return wet_wm*1.0383;
 		}
	if(Estage[2]>thr2) { if(Estage[2]>87.) return 1000;
	wet_wm=myP4(p2,Estage[2]);
	if(Estage[1]<21 || Estage[0]<19) return -1000; 
	if(Estage[2]>70.05) {
 	wet_pr=myP4(p1c,Estage[1]);
 	wet_wm=wet_wm*.3+wet_pr*.7; 
 		    }	
	return wet_wm*1.0383;
		}
	if(Estage[1]>thr1) { if(Estage[1]>87.) return 1000;
	wet_wm=myP4(p1,Estage[1]);
 	if(Estage[1]<21 || Estage[0]<19) return -1000;
	if(Estage[1]>69.7) {
	wet_pr=myP4(p0c,Estage[0]);
	wet_wm=wet_wm*.3+wet_pr*.7;	 	 
		   }    
	return wet_wm*1.0383;
		}
	if (Estage[0]>thr0) { if(Estage[1]>87.) return 1000;
	wet_wm=myP4(p0,Estage[0]);
	return wet_wm*1.0383;   // polystyrene to water equivalent
		}
	return 2000;
 }
Float_t Wepl::myP4(Float_t * p,Float_t x) { 
	Float_t x2=x*x;
return p[0]+p[1]*x+p[2]*x2+p[3]*x2*x+p[4]*x2*x2;
}
void Wepl::SetEthresholds(Float_t t0, Float_t t1, Float_t t2, Float_t t3, Float_t t4) {
   thr0=t0;thr1=t1;thr2=t2;thr3=t3;thr4=t4; 
}
