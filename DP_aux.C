/*Código desarrollado para graficar secciones eficaces y funciones usadas en el cálculo de límites de exclusión.
Última actualización: 06/04/2026
Autora: E.Depaoli
Basado en "Detecting Dark Photons with Reactor Neutrino Experiments" H.K.Park DOI: 10.1103/PhysRevLett.119.081801
*/

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <TROOT.h>
#include "TMath.h"

using std::cout;
using std::endl;
//........................................................................................
//Funciones. Explicación al final, donde están definidas.
void archivo_a_vectores(vector<double>& x_v, vector<double>& y_v, vector<double>& z_v, const string filetitle);
vector<double> gamma_ray_flux_FRJ1(vector<double> E_g_v, const double P);
vector<double> differential_cross_section_Thomson(vector<double> E_g_v, double Egprima);
vector<double> differential_cross_section_Thomson(vector<double> E_g_v, double EA, double mAc2);//DP
vector<double> differential_cross_section_Klein_Nishina(vector<double> E_g_v, double Egprima);
vector<double> differential_cross_section_Klein_Nishina(vector<double> E_g_v, double EA, double mAc2);//DP
vector<double> differential_cs_gamma_A(vector<double> dSigma_dEA_v, double epsilon);//equation (2)
vector<double> decay_length_LA(vector<double> EA_v, double mA_v, double epsilon);//(5)[mA]=[EA]=MeV
vector<double> tensor_to_vector(vector<vector<vector<double>>> LA_tensor, int ii, int jj);
//........................................................................................

//Configuracion global para los gráficos --------------------------
void SetGlobalStyle(){
	TStyle *st1 = new TStyle("st1","my style");
	//st1->SetTitleSize(0.05,"XYZ");
	//st1->SetLabelSize(0.05, "XYZ");
	st1->SetMarkerSize(1.5);
//	st1->SetPadBorderMode(0);
//	st1->SetLegendTextSize(0.035);
//	st1->SetLegendBorderSize(-1);
//	st1->SetGridColor(0);
	st1->SetLegendFillColor(0);
//	gStyle->SetLegendFillColor(0);
	gROOT->SetStyle("st1");//
	gROOT->ForceStyle();//esto no funciona: el estilo se aplica a todos los objetos creados después
	st1->cd ();//This is now the current style gStyle
}


//........................................................................................
//Contenedores···
const string filetitle{"/home/eliana/Documentos/Scripts/limites_de_exclusion/Uranium-Coef_atenuacion_nist_gov_XrayMassCoef.txt"};
vector<double> E_cross_section_v,  muroh_v, murohen_v;
vector<double> total_cross_section_v;
vector<double> E_gamma_cut_v, muroh_cut_v, murohen_cut_v;//incident gamma energy, mass atenuation coefficents. From table.
vector<double> dSgA_dE_v;
vector<double> EA_v;
vector<double> LA_v;//auxiliar. Para llenar la matriz.
vector<vector<vector<double>>> LA_tensor;
//Constantes···
float rho_U{1.895e1};//[g/cm3]
float NA{6.02214076e23};
float Mr_U{238.05079};//[g/mol]
float mec2{0.51099895000};//masa e- [MeV] 2018 CODATA value
float ro{2.8179403205e-15};//radio clásico de e- e2/4πepsilon_0mec2  [m2]
float pi_ro2mec2{0};
//........................................................................................

void DP_aux()
{
	SetGlobalStyle();
	//Coeficientes de atenuación másicos en función de la energía del gamma incidente ····
	archivo_a_vectores(E_cross_section_v,muroh_v,murohen_v,filetitle);
	
	//Corto los vectores en el intervalo de E_gamma >= 1 MeV ····
	for(int i{54}; i < E_cross_section_v.size()-1;i++)
	{
		E_gamma_cut_v.push_back(E_cross_section_v.at(i));
		muroh_cut_v.push_back(muroh_v.at(i));
		//murohen_cut_v.push_back(murohen_v.at(i));
		//cout << E_cross_section_v.at(i) << endl;
	}

	//for (auto k:E_gamma_cut_v) cout << k << endl;//mostrar en pantalla
	
	//Cuentitas ····
	pi_ro2mec2 = pow(ro,2)*TMath::Pi()*mec2;
	//cout << "pi* ro^2 * mec^2 = " << pi_ro2mec2 << endl;

	//Sección eficaz total en función de la energía del gamma incidente ···· 
	for (int i{0}; i < muroh_cut_v.size();i++) total_cross_section_v.push_back(Mr_U*muroh_cut_v.at(i)/NA/1e-24); 

	//Lleno EA_v vector de energías de los DP producidos ····
	int n{5};//cantidad de puntos en el vector
	float EA1{0.1};float EA2{4.0};//límites de integración en (7) y en la figura 1
	for(int i{0}; i < n; ++i) EA_v.push_back(i*(EA2-EA1)/n+EA1);
	//for (auto k:EA_v) cout << k << endl;//mostrar en pantalla
	cout << "#elementos en EA_v = " << EA_v.size() << endl;
	//Parámetros ····
	vector<double> epsilon_v{0.01,1.0};
	vector<double> mA_v{0.1, 0.5, 1.0};

	//····························································
	//Longitud de decaimiento de A' a 3 fotones visibles [m]  ····
	//····························································
	//En función de la energía y de la masa del A' incidente ···· 
	//Defino la dimensión del tensor vacío antes de llenarlo para poder iterar LA_tensor[i][k][j] = LA_v[k];
	//Queda LA[mA][EA][epsilon]
	vector<vector<vector<double>>> LA_tensor(mA_v.size(),vector<vector<double>>(EA_v.size(),vector<double>(epsilon_v.size())));

	//Lleno tensor 
	for(int i{0}; i < mA_v.size(); ++i){
		cout << "mA_v.at(i)" << mA_v.at(i) << " ---------- " << endl;
		for(int j{0}; j < epsilon_v.size(); ++j){
			cout << "epsilon_v.at(j)" << epsilon_v.at(j) << " ---------- " << endl;
			//Calculo LA para mA(i), epsilon(j)
			LA_v=decay_length_LA(EA_v, mA_v.at(i),epsilon_v.at(j));
	//		cout << "LA_v.size() = " << LA_v.size() << endl;
			//Guardo LA en la fila i, capa j
			for(int k{0}; k < LA_v.size();++k){
				LA_tensor[i][k][j] = LA_v[k];//
				cout << LA_v[k] << endl;
			}

			cout << "------" << endl;
			//LA_v.clear();//No es necesario porque lo piso cada vez
		}
	}

	cout << "---------------------------------------------------" << endl;
	for(int k{0}; k < EA_v.size(); ++k) cout << LA_tensor[1][k][0] << endl;
	cout << "LA_tensor size: " << std::size(LA_tensor) << endl; 
	cout << "---------------------------------------------------" << endl;

	//····························································
	//Probando Sección eficaz de producción de DP vía Compton ····
	//····························································
	//En función de la energía del gamma incidente ···· 
	vector<double> dSigma_dEA_T_v=differential_cross_section_Thomson(E_gamma_cut_v,2, 0.1);
	vector<double> dSigma_dEA_KN_v=differential_cross_section_Klein_Nishina(E_gamma_cut_v,2, 0.1);
	vector<double> dSigma_dEA_T_2_v=differential_cross_section_Thomson(E_gamma_cut_v,2, 0.5);
	vector<double> dSigma_dEA_KN_2_v=differential_cross_section_Klein_Nishina(E_gamma_cut_v,2, 0.5);
	vector<double> dSigma_dEA_T_3_v=differential_cross_section_Thomson(E_gamma_cut_v,2, 1.0);
	vector<double> dSigma_dEA_KN_3_v=differential_cross_section_Klein_Nishina(E_gamma_cut_v,2, 1.0);

	//··········································
	//Probando sección eficaz Compton usual ····
	//··········································
	//En función de la energía del gamma incidente ···· 
	vector<double> dSigma_dEgprima_T_v=differential_cross_section_Thomson(E_gamma_cut_v,5.11e-3);
	vector<double> dSigma_dEgprima_KN_v=differential_cross_section_Klein_Nishina(E_gamma_cut_v,5.11e-3);

	vector<double> dSigma_dEgprima_T_2_v=differential_cross_section_Thomson(E_gamma_cut_v,5.11e-2);
	vector<double> dSigma_dEgprima_KN_2_v=differential_cross_section_Klein_Nishina(E_gamma_cut_v,5.11e-2);

	vector<double> dSigma_dEgprima_T_3_v=differential_cross_section_Thomson(E_gamma_cut_v,5.11e-1);
	vector<double> dSigma_dEgprima_KN_3_v=differential_cross_section_Klein_Nishina(E_gamma_cut_v,5.11e-1);
	//···································

	//Flujo gamma generado en el reactor en función de la energía ···· 
	vector<double> dNg_dEg_v=gamma_ray_flux_FRJ1(E_gamma_cut_v,1000);

	//Flujo gamma total generado en el reactor ···· 
	float integral_above_1MeV{0};
	
	for(int i{0}; i < E_gamma_cut_v.size()-1;i++)
	{
		integral_above_1MeV += dNg_dEg_v.at(i)*(E_gamma_cut_v.at(i+1)-E_gamma_cut_v.at(i));
		//cout << "dE_g = " << E_gamma_cut_v.at(i+1) -E_gamma_cut_v.at(i) << " -- Integral acumulada = " << integral_above_1MeV << endl;//
	}

	//····························································
	//Producción de DP vía Compton 							  ····
	//····························································
	vector<double> dSgamma_A_v=differential_cs_gamma_A(dSigma_dEA_T_v, 1);///(E_gamma_cut_v,2, 0.1);//differential cross section for sigma gamma -> A'

	/*
	cout << "Longitud(energías) = " << E_gamma_cut_v.size() << endl;
	cout << "Longitud(flujo gamma) = " << dNg_dEg_v.size() << endl;
	cout << "Longitud(sección eficaz total) = " << total_cross_section_v.size() << endl;
	cout << "Longitud(sección eficaz diferencial gamma -> A') = " << dSgamma_A_v.size() << endl;
*/

	//Calculo un punto en la figura 1. dNA'/dEA'
/*	double dNA_dEA{0};

	for(int i{0}; i < dSgamma_A_v.size()-1; ++i){
		dNA_dEA+= pow(total_cross_section_v.at(i),-1)*dSgamma_A_v.at(i)*dNg_dEg_v.at(i)*(E_gamma_cut_v.at(i+1)-E_gamma_cut_v.at(i));
	}
*/
	//for (auto k:dSgamma_A_v) cout << k << endl;//mostrar en pantalla
	
	//Gráficos ····
	//····························································
	//Longitud de decaimiento de A' a 3 fotones visibles [m]  ····
	//····························································
	//Vectores auxiliares para graficar mi "tensor"
	vector<double> La_v1=tensor_to_vector(LA_tensor, 0, 1);
	vector<double> La_v2=tensor_to_vector(LA_tensor, 1, 1);
	vector<double> La_v3=tensor_to_vector(LA_tensor, 2, 1);
	vector<double> La_v11=tensor_to_vector(LA_tensor, 0, 0);
	vector<double> La_v22=tensor_to_vector(LA_tensor, 1, 0);
	vector<double> La_v33=tensor_to_vector(LA_tensor, 2, 0);
	
	TCanvas*c5 = new TCanvas("DP decay length","DP decay length", 900, 500);
	TGraph*LA_v = new TGraph(EA_v.size(), &EA_v[0], &La_v1[0]);
	TGraph*LA11_v = new TGraph(EA_v.size(), &EA_v[0], &La_v11[0]);
	TGraph*LA2_v = new TGraph(EA_v.size(), &EA_v[0], &La_v2[0]);
	TGraph*LA22_v = new TGraph(EA_v.size(), &EA_v[0], &La_v22[0]);
	TGraph*LA3_v = new TGraph(EA_v.size(), &EA_v[0], &La_v3[0]);
	TGraph*LA33_v = new TGraph(EA_v.size(), &EA_v[0], &La_v33[0]);
	
	c5->SetGrid();
	gPad->SetLogy();
	gPad->SetLogx();
	LA_v->Draw("AP");
	LA2_v->Draw("sameP");//
	LA3_v->Draw("sameP");//
	LA11_v->Draw("sameP");//
	LA22_v->Draw("sameP");//
	LA33_v->Draw("sameP");//
	
	LA_v->SetMarkerStyle(20);
	LA2_v->SetMarkerStyle(22);
	LA3_v->SetMarkerStyle(21);
	LA11_v->SetMarkerStyle(20);
	LA22_v->SetMarkerStyle(22);
	LA33_v->SetMarkerStyle(21);

	LA_v->SetMarkerColor(kCyan);
	LA2_v->SetMarkerColor(kCyan+1);
	LA3_v->SetMarkerColor(kBlue-4);
	
	LA11_v->SetMarkerColor(kOrange);
	LA22_v->SetMarkerColor(kOrange+1);
	LA33_v->SetMarkerColor(kOrange-4);
	
	LA_v->SetTitle("Dark Photon decay length");
	LA_v->GetYaxis()->SetTitle("L_{A'}[m]");
	LA_v->GetYaxis()->SetRangeUser(1e1, 1e18);
	LA_v->GetXaxis()->SetTitle("E [MeV]");

	auto leg5 = new TLegend(0.70,0.67,0.90,0.97); 
    leg5->AddEntry(LA_v,("m'_{A'} = "+std::to_string(mA_v.at(0)).substr(0,3)+" MeV " + "#varepsilon = " + std::to_string(epsilon_v.at(1)).substr(0,4)).c_str(),"p");
    leg5->AddEntry(LA2_v,("m'_{A'} = "+std::to_string(mA_v.at(1)).substr(0,3)+" MeV "+ "#varepsilon = " + std::to_string(epsilon_v.at(1)).substr(0,4)).c_str(),"p");
    leg5->AddEntry(LA3_v,("m'_{A'} = "+std::to_string(mA_v.at(2)).substr(0,3)+" MeV "+ "#varepsilon = " + std::to_string(epsilon_v.at(1)).substr(0,4)).c_str(),"p");
    leg5->AddEntry(LA11_v,("m'_{A'} = "+std::to_string(mA_v.at(0)).substr(0,3)+" MeV " + "#varepsilon = " + std::to_string(epsilon_v.at(0)).substr(0,4)).c_str(),"p");
    leg5->AddEntry(LA22_v,("m'_{A'} = "+std::to_string(mA_v.at(1)).substr(0,3)+" MeV "+ "#varepsilon = " + std::to_string(epsilon_v.at(0)).substr(0,4)).c_str(),"p");
    leg5->AddEntry(LA33_v,("m'_{A'} = "+std::to_string(mA_v.at(2)).substr(0,3)+" MeV "+ "#varepsilon = " + std::to_string(epsilon_v.at(0)).substr(0,4)).c_str(),"p");

    leg5->SetTextSize(0.03); 
    leg5->SetBorderSize(0);
    leg5->Draw();

    //c5->Print("DP_decay_length.png");
	//------------------------------------------------------------

	//------------------------------------------------------------
	//Secciones eficaces de producción
	//------------------------------------------------------------
/*	TCanvas*c4 = new TCanvas("DP production differential cross section","DP production differential cross section", 900, 500);
	TGraph*TDP_gr = new TGraph(E_gamma_cut_v.size(), &E_gamma_cut_v[0], &dSigma_dEA_T_v[0]);
	TGraph*KNDP_gr = new TGraph(E_gamma_cut_v.size(), &E_gamma_cut_v[0], &dSigma_dEA_KN_v[0]);

	TGraph*TDP2_gr = new TGraph(E_gamma_cut_v.size(), &E_gamma_cut_v[0], &dSigma_dEA_T_2_v[0]);
	TGraph*KNDP2_gr = new TGraph(E_gamma_cut_v.size(), &E_gamma_cut_v[0], &dSigma_dEA_KN_2_v[0]);
	
	TGraph*TDP3_gr = new TGraph(E_gamma_cut_v.size(), &E_gamma_cut_v[0], &dSigma_dEA_T_3_v[0]);
	TGraph*KNDP3_gr = new TGraph(E_gamma_cut_v.size(), &E_gamma_cut_v[0], &dSigma_dEA_KN_3_v[0]);
	
	c4->SetGrid();
	gPad->SetLogy();
	gPad->SetLogx();
	//TDP_gr->Draw("AP");
	//KNDP_gr->Draw("sameP");//
	TDP_gr->SetMarkerStyle(20);
	KNDP_gr->SetMarkerStyle(20);

	TDP2_gr->Draw("sameP");//
	KNDP2_gr->Draw("sameP");//
	TDP3_gr->Draw("sameP");//
	KNDP3_gr->Draw("sameP");//

	TDP2_gr->SetMarkerStyle(22);
	KNDP2_gr->SetMarkerStyle(22);
	TDP3_gr->SetMarkerStyle(21);
	KNDP3_gr->SetMarkerStyle(21);
	
	TDP_gr->SetMarkerColor(kCyan);
	KNDP_gr->SetMarkerColor(kOrange);

	TDP2_gr->SetMarkerColor(kCyan+1);
	KNDP2_gr->SetMarkerColor(kOrange+7);
	TDP3_gr->SetMarkerColor(kBlue-4);
	KNDP3_gr->SetMarkerColor(kOrange+10);
	
	TDP_gr->SetTitle("Differential cross section E_{A'} = 2.0 MeV ");
	TDP_gr->GetYaxis()->SetTitle("#frac{d#sigma_{C}}{dE_{A'}}(E'_{A'})[m^{2}MeV^{-1}s^{-1}]");
	TDP_gr->GetYaxis()->SetRangeUser(1e-32, 1e-28);
	TDP_gr->GetXaxis()->SetTitle("E [MeV]");

	auto leg4 = new TLegend(0.70,0.67,0.90,0.97); 
    leg4->AddEntry(TDP_gr,("Thomson m'_{A'} = "+std::to_string(0.1).substr(0,3)+" MeV").c_str(),"p");
    leg4->AddEntry(KNDP_gr,("Klein-Nishina m'_{A'} = "+std::to_string(0.1).substr(0,3)+" MeV").c_str(),"p");
    leg4->AddEntry(TDP2_gr,("Thomson m'_{A'} = "+std::to_string(0.5).substr(0,3)+" MeV").c_str(),"p");
    leg4->AddEntry(KNDP2_gr,("Klein-Nishina m'_{A'} = "+std::to_string(0.5).substr(0,3)+" MeV").c_str(),"p");
    leg4->AddEntry(TDP3_gr,("Thomson m'_{A'} = "+std::to_string(1.0).substr(0,3)+" MeV").c_str(),"p");
    leg4->AddEntry(KNDP3_gr,("Klein-Nishina m'_{A'} = "+std::to_string(1.0).substr(0,3)+" MeV").c_str(),"p");

    leg4->SetTextSize(0.03); 
    leg4->SetBorderSize(0);
    leg4->Draw();

    //c4->Print("DP_production_differential_cross_section.png");
	//------------------------------------------------------------
	//Secciones eficaces de producción
	//------------------------------------------------------------
	TCanvas*c3 = new TCanvas("Usual Compton differential cross section","Usual Compton differential cross section", 900, 500);
	TGraph*Thomson_gr = new TGraph(E_gamma_cut_v.size(), &E_gamma_cut_v[0], &dSigma_dEgprima_T_v[0]);
	TGraph*KN_gr = new TGraph(E_gamma_cut_v.size(), &E_gamma_cut_v[0], &dSigma_dEgprima_KN_v[0]);

	TGraph*T2_gr = new TGraph(E_gamma_cut_v.size(), &E_gamma_cut_v[0], &dSigma_dEgprima_T_2_v[0]);
	TGraph*KN2_gr = new TGraph(E_gamma_cut_v.size(), &E_gamma_cut_v[0], &dSigma_dEgprima_KN_2_v[0]);
	
	TGraph*T3_gr = new TGraph(E_gamma_cut_v.size(), &E_gamma_cut_v[0], &dSigma_dEgprima_T_3_v[0]);
	TGraph*KN3_gr = new TGraph(E_gamma_cut_v.size(), &E_gamma_cut_v[0], &dSigma_dEgprima_KN_3_v[0]);
	
	c3->SetGrid();
	gPad->SetLogy();
	gPad->SetLogx();
	Thomson_gr->Draw("AP");
	KN_gr->Draw("sameP");//

	T2_gr->Draw("sameP");//
	KN2_gr->Draw("sameP");//
	T3_gr->Draw("sameP");//
	KN3_gr->Draw("sameP");//

	Thomson_gr->SetMarkerStyle(20);
	KN_gr->SetMarkerStyle(20);
	
	T2_gr->SetMarkerStyle(22);
	KN2_gr->SetMarkerStyle(22);
	T3_gr->SetMarkerStyle(21);
	KN3_gr->SetMarkerStyle(21);
	
	Thomson_gr->SetMarkerColor(kCyan);
	KN_gr->SetMarkerColor(kOrange);

	T2_gr->SetMarkerColor(kCyan+1);
	KN2_gr->SetMarkerColor(kOrange+7);
	T3_gr->SetMarkerColor(kBlue-4);
	KN3_gr->SetMarkerColor(kOrange+10);
	
	Thomson_gr->SetTitle("Differential cross section");
	Thomson_gr->GetYaxis()->SetTitle("#frac{d#sigma_{C}}{dE_{#gamma}}(E'_{#gamma})[m^{2}MeV^{-1}s^{-1}]");
	Thomson_gr->GetYaxis()->SetRangeUser(1e-32, 1e-14);
	Thomson_gr->GetXaxis()->SetTitle("E [MeV]");

	auto leg3 = new TLegend(0.70,0.67,0.90,0.97); 
    leg3->AddEntry(Thomson_gr,("Thomson E'_{#gamma} = "+std::to_string(5.11).substr(0,4)+" keV").c_str(),"p");
    leg3->AddEntry(KN_gr,("Klein-Nishina E'_{#gamma} = "+std::to_string(5.11).substr(0,4)+" keV").c_str(),"p");
    leg3->AddEntry(T2_gr,("Thomson E'_{#gamma} = "+std::to_string(5.11e1).substr(0,4)+" keV").c_str(),"p");
    leg3->AddEntry(KN2_gr,("Klein-Nishina E'_{#gamma} = "+std::to_string(5.11e1).substr(0,4)+" keV").c_str(),"p");
    leg3->AddEntry(T3_gr,("Thomson E'_{#gamma} = "+std::to_string(5.11e2).substr(0,3)+" keV").c_str(),"p");
    leg3->AddEntry(KN3_gr,("Klein-Nishina E'_{#gamma} = "+std::to_string(5.11e2).substr(0,3)+" keV").c_str(),"p");

    leg3->SetTextSize(0.03); 
    leg3->SetBorderSize(0);
    leg3->Draw();
    //c3->Print("Usual_Compton_differential_cross_section.png");
	//------------------------------------------------------------
	//Producción de gammas
	//------------------------------------------------------------

	TCanvas*c2 = new TCanvas("Gamma flux FRJ-1","Gamma flux FRJ-1", 900, 500);
	TGraph*gfFRJ = new TGraph(E_gamma_cut_v.size(), &E_gamma_cut_v[0], &dNg_dEg_v[0]);
	gfFRJ->Draw("AP");
	gfFRJ->SetMarkerStyle(20);
	gfFRJ->SetMarkerColor(kBlue);
	c2->SetGrid();
	gPad->SetLogy();
	gPad->SetLogx();
	gfFRJ->SetTitle("Prompt gamma ray");
	gfFRJ->GetYaxis()->SetTitle("#frac{dN_{#gamma}}{dE_{#gamma}}[MeV^{-1}s^{-1}]");
	gfFRJ->GetYaxis()->SetRangeUser(1e11, 1e21);
	gfFRJ->GetXaxis()->SetTitle("E [MeV]");

	//c2->Print("FRJ-1reactor_gamma_flux.png");
	//------------------------------------------------------------
	//Sección eficaz total de interacción de gammas contra Uranio
	//------------------------------------------------------------
	TCanvas*c1 = new TCanvas("Mass Attenuation Coeficients","Mass Attenuation Coeficients", 900, 500);
	TGraph*mac = new TGraph(E_gamma_cut_v.size(), &E_gamma_cut_v[0], &muroh_cut_v[0]);//&E_gamma_cut_v[0]);//
	TGraph*csU = new TGraph(E_gamma_cut_v.size(), &E_gamma_cut_v[0], &total_cross_section_v[0]);
	c1->SetGrid();
	gPad->SetLogy();
	gPad->SetLogx();
	//mac->GetXaxis()->SetTitle("E [MeV]");
	csU->GetXaxis()->SetTitle("E [MeV]");
	//mac->GetYaxis()->SetTitle("#mu/#rho [cm2/g]");
	//csU->GetYaxis()->SetTitle("[cm2]");
	csU->GetYaxis()->SetTitle("#sigma [barn]");
	csU->SetTitle("Uranium");
	mac->SetMarkerStyle(20);
	csU->SetMarkerStyle(22);
	mac->SetMarkerColor(kBlue);
	csU->SetMarkerColor(kOrange+2);
	//mac->GetYaxis()->SetRangeUser(1e-2, 1e27);
	//mac->Draw("AP");
	csU->Draw("AP");//Draw("same");//
	//c1->Print("total_cross_section.png");
*/
	// Fin void DP ····
}

void archivo_a_vectores(vector<double>& x_v, vector<double>& y_v, vector<double>& z_v, const string filetitle)
{
	
	std::ifstream file;//guardará temporalmente el contenido del archivo
	file.open(filetitle);

	//En caso de que no pueda abrir o localizar el archivo
	if(!file.is_open())
	{
		std::cerr << filetitle << "No puedo abrirlo\n";
		return;
	}

	std::string stream;
	for(int i=0; i < 11; i++) std::getline(file,stream);

	double x,y,z;
	while(!file.eof())
	{
		file >> x >> y >> z;//
		x_v.push_back(x);
		y_v.push_back(y);
		z_v.push_back(z);
	}

	file.close();
}

//Guardar las filas del tensor [i][k][j] en un vector para graficar
vector<double> tensor_to_vector(vector<vector<vector<double>>> LA_tensor, int ii, int jj){
	vector<double> y_val;
	for(int k = 0; k < EA_v.size(); ++k) {
	    y_val.push_back(LA_tensor[ii][k][jj]);
	}
	return y_val;
}

//DP decay length. Eq. (5)
//Use [EA]=MeV, [mA]=MeV
vector<double> decay_length_LA(vector<double> EA_v, double mA_v, double epsilon){
	vector<double> LA_v;
	for(auto k:EA_v) LA_v.push_back(505*pow(epsilon,-2)*k*pow(mA_v,-10));
	return LA_v;
}

//Gamma flux for the FRJ-1 reactor, valid for E_gamma > 200 keV. Use [Energy]=MeV [P] = MW
vector<double> gamma_ray_flux_FRJ1(vector<double> E_g_v, const double P){
	vector<double> dNg_dEg_v;
	float aux{0.91};//si este se achica, el valor de energía donde se produce el codo se mueve hacia la izquierda. Para 0.91, el codo representa una caída de 63% del valor inicial en 1 MeV
	for(auto k: E_g_v) dNg_dEg_v.push_back(0.58e18*P*TMath::Exp(-k/aux));
	return dNg_dEg_v;	
}

//differential cross section for sigma gamma -> A' in the limit mA' << me
vector<double> differential_cs_gamma_A(vector<double> dSigma_dEA_v, double epsilon){
	vector<double> dSgA_dEA_v;
	for(auto k:dSigma_dEA_v) dSgA_dEA_v.push_back(pow(epsilon,2)*k);
	return dSgA_dEA_v;
}

//Interacción tipo Compton: gamma + e- -> A + e-
//Sección eficaz diferencial en función de la energía del gamma incidente para un valor fijo de energía y masa del DP dispersado
//Use [Energy]=MeV
//Thomson····
vector<double> differential_cross_section_Thomson(vector<double> E_g_v, double EA, double mAc2){
	vector<double> dSigma_dEA_v;
	double num{0};
	double denom{0};
	double cos_2_theta{0};
	
	for(auto k: E_g_v){
		num=pow(k*mec2+pow(mAc2,2)/2-(k+mec2)*EA,2);
		denom=pow(k,2)*(pow(EA,2)-pow(mAc2,2));
		cos_2_theta=num/denom;
		dSigma_dEA_v.push_back((cos_2_theta+1)*pow(EA,-2)*pi_ro2mec2);
		//cout << "cos_2_theta = "<< cos_2_theta << endl;//"num = " << num << "denom = " << denom <<
		//cout << (cos_2_theta+1)*pow(EA,-2)*pi_ro2mec2 << endl;
	}

	return dSigma_dEA_v;
}
//Klein Nishina····
vector<double> differential_cross_section_Klein_Nishina(vector<double> E_g_v, double EA, double mAc2){
	vector<double> dSigma_dEA_v;
	double num{0};
	double denom{0};
	double cos_2_theta{0};

	for(auto k: E_g_v){
		num=pow(k*mec2+pow(mAc2,2)/2-(k+mec2)*EA,2);
		denom=pow(k,2)*(pow(EA,2)-pow(mAc2,2));
		cos_2_theta=num/denom;
		dSigma_dEA_v.push_back((cos_2_theta + k*pow(EA,-1) + EA*pow(k,-1) - 1)*pow(k,-3)*EA*pi_ro2mec2);
		//cout << "cos_2_theta = "<< cos_2_theta << endl;
	}

	return dSigma_dEA_v;
}

//Compton clásico
//Sección eficaz diferencial en función de la energía del gamma incidente para un valor fijo de energía del gamma dispersado
//Use [Energy]=MeV
//Thomson····
vector<double> differential_cross_section_Thomson(vector<double> E_g_v, double Egprima){
	vector<double> dSigma_dEgprima_v;
	double cos_2_theta{0};
	for(auto k: E_g_v){
		cos_2_theta=pow(mec2*(pow(k,-1)-pow(Egprima, -1))+1,2);
		dSigma_dEgprima_v.push_back((cos_2_theta+1)*pow(Egprima,-2)*pi_ro2mec2);
	}
	return dSigma_dEgprima_v;
}
//Klein Nishina····
vector<double> differential_cross_section_Klein_Nishina(vector<double> E_g_v, double Egprima){
	vector<double> dSigma_dEgprima_v;
	double cos_2_theta{0};
	for(auto k: E_g_v){
		cos_2_theta=pow(mec2*(pow(k,-1)-pow(Egprima, -1))+1,2);
		dSigma_dEgprima_v.push_back((cos_2_theta + k*pow(Egprima,-1) + Egprima*pow(k,-1) - 1)*pow(k,-3)*Egprima*pi_ro2mec2);
	}
	return dSigma_dEgprima_v;
}

//cout << "Hola" << endl;