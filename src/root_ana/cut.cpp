#include <iostream>
#include <stdlib.h> 

#include "TString.h"


// /*--------------------------------------------------*/
// Define Cuts
// Define HMS Cuts
//
/*--------------------------------------------------*/
// Kumac code HMS cut from Paw script
//    cuts $11 abs(hsdelta)<8.0
// *   cuts $12 abs(hsxptar)<0.060  | cut on Jochen's thesis p.53
//    cuts $12 abs(hsxptar)<0.080   | extra sieve slit fitting should improve things
//    cuts $13 abs(hsyptar)<0.035
// 
//    cuts $10 $11.and.$12.and.$13

using namespace std;


TString HMS_cut() {

	TString hms_cut = "";

 	return hms_cut;

}

/*--------------------------------------------------*/
// /*--------------------------------------------------*/
// Define SOS Cuts from PAW script
//    cuts $21 abs(ssdelta)<15                       
//    cuts $22 abs(ssxfp)<20	
//    cuts $23 abs(ssxfp+ssxpfp*313.01+5.64)<50.
//    cuts $24 abs(ssyfp+ssypfp*313.01)<30.
//    cuts $25 ssytar<1.5          | cut on Jochen's thesis p.80, for th_SOS>42deg.
//    cuts $26 abs(ssxptar)<0.04
//    cuts $27 abs(ssyptar)<0.065                 

//    Jochen acceptance cut (see p.17 of Blok paper)
//    t1=-0.125+0.00425*ssdelta+0.064*ssytar-0.0017*ssdelta*ssytar
//    t2= 0.125-0.00425*ssdelta+0.064*ssytar-0.0017*ssdelta*ssytar 
//    cuts $28 [t1]<ssyptar<[t2]
//
//    cuts $20 $21.and.$22.and.$23.and.$24.and.$25.and.$26.and.$27.and.$28

TString SOS_cut() {

	TString sos_cut = "";

   	return sos_cut; 

}



/*--------------------------------------------------*/
// /*--------------------------------------------------*/
// Define PID Cuts from PAW Script
// 

//    cuts $31 hsbeta>0.95
//    cuts $32 (haero_po+haero_ne)>3.0
//    cuts $33 hcer_npe<2.
//    cuts $34 scer_npe>0.5
//    cuts $35 ssshtrk>0.70           
// 
//    cuts $30 $31.and.$32.and.$33.and.$34.and.$35



TString PID_cut() {

	TString pid_cut = "";

	return pid_cut;

}





/*--------------------------------------------------*/
// /*--------------------------------------------------*/
// Kumac diamand cuts from PAW Script
//
//    cuts $51 q2>(-2.916*(w-2.4733)+1.5846)
//    cuts $52 q2>(-4.452*(w-1.94)+3.4199) 
//    cuts $53 q2<(-2.677*(w-2.3308)+2.231)
//    cuts $54 q2<(-4.3748*(w-2.41)+1.78)


TString Diamond_cut() {

 	TString cut1 = ""; 

	return cut1;

}





TString Missingmass_cut(Double_t m_m_offset) {

	TString missmass_cut = ""; 

	return missmass_cut;

}


TString Set_t_limit(Double_t tmin, Double_t tmax) {

	TString t_limit = ""; 

	return t_limit;
}




TString Cointime_primary_cut(Double_t center) {

	TString cut_tmp= "";

	return cut_tmp;
}




TString Cointime_random_cut(Double_t center) {

 	TString cut1 = ""; 

	return cut1;
}


TString Cointime_all(Double_t center) {


 	TString cut1 = ""; 

	return cut1;
}







TString Heep_PID_Cut(float Q2_set) {

	TString pid_cut = "";

	return pid_cut;

}


TString T_heep_cut() {

	TString missmass_cut = ""; 

	return missmass_cut;

}





TString Missingmass_heep_cut() {

	TString total_cut = "";

	return total_cut;

}




/*--------------------------------------------------*/
/*--------------------------------------------------*/
/// Singles check cut on W

TString Heep_single_W_cut(float w_down, float w_up) {

	TString missmass_cut = ""; 

	return missmass_cut;

}

/*--------------------------------------------------*/
/*--------------------------------------------------*/
/// Singles PID cut

TString Heep_Single_PID_cut() {

	TString pid_cut = "";

	return pid_cut;

} 


/*--------------------------------------------------*/
/*--------------------------------------------------*/
/// Singles SOS cut
TString Heep_Single_SOS_cut() {

	TString sos_cut = "";

   	return sos_cut; 



}


/*--------------------------------------------------*/
/*--------------------------------------------------*/
/// Garth's Omega cuts form teleconf report dated 05/April/2005

TString Omega_HMS_Cut_Garth() {

	TString hms_cut = "";

 	return hms_cut;

}

TString Omega_SOS_Cut_Garth() {

	TString sos_cut = "";

   	return sos_cut; 

}

TString Omega_PID_Cut_Garth() {

	TString pid_cut = "";

	return pid_cut;

}


/*--------------------------------------------------*/
/*--------------------------------------------------*/
/// Modified cut based on Garth's Omega cuts and stardard pion analysis cut
//  01/March/2016

TString Omega_HMS_Cut() {


	TString hms_cut = "";

 	return hms_cut;

}


TString Omega_SOS_Cut() {


	TString sos_cut = "";

   	return sos_cut; 

}


TString Omega_PID_Cut() {

	TString pid_cut = "";

	return pid_cut;

}


TString Omega_U_Cut() {

 	TString u_cut = "";
 
	return u_cut;

}



TString Omega_U_Sim_Cut() {

 	TString u_cut = "";
 
	return u_cut;

}



/*--------------------------------------------------*/
// Missmass cut

TString Omega_Missmass_Cut() {

 	TString mm_cut = "";
 
	return mm_cut;

}

