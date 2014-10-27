// Class CalF members are constructor CalF and  function CalTVf 
// converting  ADC values into energy deposited in a stage (MeV)
// with 2D spatial correction utilizing T and V coordinates 
// of the energy detector hit indicated by track extrapolation.
//Created by V.A.Bashkirov, May 9,2014
 using std::cout;     
 using std::endl;
 class CalF {  
 Float_t p0,p1,p2,p3,p4;    // po-p4 store T-V fit parameters
 Float_t pn;  // calibration coefficient equal to average energy  
	      // deposited in corresponding stage by 200MeV proton  
  	      // to convert ADC value into energy deposition (MeV)
  public: 
	CalF(Int_t, Float_t* );  //constructor 
	Float_t CalTVf(Float_t,Float_t,Float_t); //function
};
CalF::CalF(Int_t Nstage, Float_t* Fpar ) {
  static Float_t Estage[5]={25.85,29.04,34.43,47.8,52.63};	   
  p0=Fpar[0]; // Function initialization with corresponding
  p1=Fpar[1]; // T-V fit parameters and calibration coefficient 
  p2=Fpar[2]; // assignment for stage # Nstage
  p3=Fpar[3]; 
  p4=Fpar[4];
  if( Nstage<5 && Nstage>=0) pn=Estage[Nstage];
  else cout<<"Wrong stage # "<<Nstage<<"! Valid range 0-4"<<endl;
}
Float_t CalF::CalTVf(Float_t t,Float_t v, Float_t adc) {
return pn*adc/(p0+p1*t+p2*t*t+p3*v+p4*v*v); // function itself
}
