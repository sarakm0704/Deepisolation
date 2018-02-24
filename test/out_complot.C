void out_complot(string inputfile){

  TFile fileIn(Form("%s",inputfile.c_str()));
  TTreeReader reader("tree",&fileIn);
  TTreeReaderArray<Double_t> chN(reader, "chN");
  TTreeReaderArray<Double_t> ch_dR(reader, "ch_dR");
  TTreeReaderArray<Double_t> ch_dEta(reader, "ch_deta");
  TTreeReaderArray<Double_t> ch_dPhi(reader, "ch_dphi");
  TTreeReaderArray<Double_t> chIso(reader, "chIso");

  TTreeReaderArray<Double_t> nhN(reader, "nhN");
  TTreeReaderArray<Double_t> nh_dR(reader, "nh_dR");
  TTreeReaderArray<Double_t> nh_dEta(reader, "nh_deta");
  TTreeReaderArray<Double_t> nh_dPhi(reader, "nh_dphi");
  TTreeReaderArray<Double_t> nhIso(reader, "nhIso");

  TTreeReaderArray<Double_t> phN(reader, "phN");
  TTreeReaderArray<Double_t> ph_dR(reader, "ph_dR");
  TTreeReaderArray<Double_t> ph_dEta(reader, "ph_deta");
  TTreeReaderArray<Double_t> ph_dPhi(reader, "ph_dphi");
  TTreeReaderArray<Double_t> phIso(reader, "phIso");
 
  TTreeReaderArray<Double_t> pT(reader, "pt");
  TTreeReaderArray<Double_t> Eta(reader, "eta");
  TTreeReaderArray<Double_t> Phi(reader, "phi");
 
  TH1D *h_chN = new TH1D("h_chN","h_chN",30,0,30);
  TH1D *h_ch_dR = new TH1D("h_ch_dR","h_ch_dR",100,0,0.5);
  TH1D *h_ch_dEta = new TH1D("h_ch_deta","h_ch_deta",100,0,0.5);
  TH1D *h_ch_dPhi = new TH1D("h_ch_dphi","h_ch_dphi",100,0,0.5);
  TH1D *h_chIso = new TH1D("h_chIso","h_chIso",100,0,50);
  
  TH1D *h_nhN = new TH1D("h_nhN","h_nhN",30,0,30);
  TH1D *h_nh_dR = new TH1D("h_nh_dR","h_nh_dR",100,0,0.5);
  TH1D *h_nh_dEta = new TH1D("h_nh_deta","h_nh_deta",100,0,0.5);
  TH1D *h_nh_dPhi = new TH1D("h_nh_dphi","h_nh_dphi",100,0,0.5);
  TH1D *h_nhIso = new TH1D("h_nhIso","h_nhIso",100,0,50);
  
  TH1D *h_phN = new TH1D("h_phN","h_phN",30,0,30);
  TH1D *h_ph_dR = new TH1D("h_ph_dR","h_ph_dR",100,0,0.5);
  TH1D *h_ph_dEta = new TH1D("h_ph_deta","h_ph_deta",100,0,0.5);
  TH1D *h_ph_dPhi = new TH1D("h_ph_dphi","h_ph_dphi",100,0,0.5);
  TH1D *h_phIso = new TH1D("h_phIso","h_phIso",100,0,50);

  TH1D *h_pT = new TH1D("h_pt","h_pt",100,0,100);
  TH1D *h_Eta = new TH1D("h_eta","h_eta",100,-6,6);
  TH1D *h_Phi = new TH1D("h_phi","h_phi",100,-5,5);

  ofstream myfile;
  myfile.open ("qcd_iso.txt");
 
  while(reader.Next()){
  
    if( chN.GetSize() == 0 ) continue;
    double chn = -1;
    //cout << "size = " << chN.GetSize() << endl;
    for (int i = 0; i < chN.GetSize(); ++i){
      chn = chN[i];
      break;
    } 
    h_chN->Fill(chn);
    myfile << chn << " ";

    double ch_dr = -1;
    double ch_deta = -1;
    double ch_dphi = -1;
    double chiso = -1;

    if( chn == 0) {
     ch_dr = 0;
     ch_deta = 0;
     ch_dphi = 0;
     chiso = 0; 
    }else{ 
      for (int i = 0; i < ch_dR.GetSize(); ++i){
        ch_dr = ch_dR[i];
        break;
      } 
      for (int i = 0; i < ch_dEta.GetSize(); ++i){
        ch_deta = ch_dEta[i];
        break;
      } 
      for (int i = 0; i < ch_dPhi.GetSize(); ++i){
        ch_dphi = ch_dPhi[i];
        break;
      }
      for (int i = 0; i < chIso.GetSize(); ++i){
        chiso = chIso[i];
        break;
      }
    }
 
    h_ch_dR->Fill(ch_dr);
    myfile << ch_dr << " ";
    h_ch_dEta->Fill(ch_deta);
    myfile << ch_deta << " "; 
    h_ch_dPhi->Fill(ch_dphi);
    myfile << ch_dphi << " ";
    h_chIso->Fill(chiso);
    myfile << chiso << " ";

    if( nhN.GetSize() == 0 ) continue;
    double nhn = -1;
    for (int i = 0; i < nhN.GetSize(); ++i){
      nhn = nhN[i];
      break;
    }
    h_nhN->Fill(nhn);
    myfile << nhn << " "; 
  
    double nh_dr = -1;
    double nh_deta = -1;
    double nh_dphi = -1;
    double nhiso = -1;

    if( nhn == 0) {
     nh_dr = 0;
     nh_deta = 0;
     nh_dphi = 0; 
     nhiso = 0;
    }else{ 
      for (int i = 0; i < nh_dR.GetSize(); ++i){
        nh_dr = nh_dR[i];
        break;
      }
      for (int i = 0; i < nh_dEta.GetSize(); ++i){
        nh_deta = nh_dEta[i];
        break;
      }
      for (int i = 0; i < nh_dPhi.GetSize(); ++i){
        nh_dphi = nh_dPhi[i];
        break;
      }
      for (int i = 0; i < nhIso.GetSize(); ++i){
        nhiso = nhIso[i];
        break;
      }
    }
    h_nh_dR->Fill(nh_dr);
    myfile << nh_dr << " ";
    h_nh_dEta->Fill(nh_deta);
    myfile << nh_deta << " ";
    h_nh_dPhi->Fill(nh_dphi);
    myfile << nh_dphi << " ";
    h_nhIso->Fill(nhiso);
    myfile << nhiso << " ";

    if( phN.GetSize() == 0 ) continue;
    double phn = -1;
    for (int i = 0; i < phN.GetSize(); ++i){
      phn = phN[i];
      break;
    }
    h_phN->Fill(phn);
    myfile << phn << " ";

    double ph_dr = -1;
    double ph_deta = -1;
    double ph_dphi = -1;
    double phiso = -1;
  
    if( phn == 0) {
     ph_dr = 0;
     ph_deta = 0;
     ph_dphi = 0; 
     phiso = 0;

    }else{ 
      for (int i = 0; i < ph_dR.GetSize(); ++i){
        ph_dr = ph_dR[i];
        break;
      }
      for (int i = 0; i < ph_dEta.GetSize(); ++i){
        ph_deta = ph_dEta[i];
        break;
      }
      for (int i = 0; i < ph_dPhi.GetSize(); ++i){
        ph_dphi = ph_dPhi[i];
        break;
      }
      for (int i = 0; i < phIso.GetSize(); ++i){
        phiso = phIso[i];
        break;
      }
    }
    h_ph_dR->Fill(ph_dr);
    myfile << ph_dr << " ";
    h_ph_dEta->Fill(ph_deta);
    myfile << ph_deta << " ";
    h_ph_dPhi->Fill(ph_dphi);
    myfile << ph_dphi << " ";
    h_phIso->Fill(phiso);
    myfile << phiso << " ";

  //  cout << "pT" << endl;

    if( pT.GetSize() == 0 ) continue;
    double pt = -1;
    for (int i = 0; i < pT.GetSize(); ++i){
      pt = pT[i];
   //   cout << pt << endl;
      break;
    }
    h_pT->Fill(pt);
    //myfile << pt << " ";
 
    if( Eta.GetSize() == 0 ) continue;
    double eta = -1;
    for (int i = 0; i < Eta.GetSize(); ++i){
      eta = Eta[i];
      break;
    }
    h_Eta->Fill(eta);
    //myfile << eta << " ";

    if( Phi.GetSize() == 0 ) continue;
    double phi = -1;
    for (int i = 0; i < Phi.GetSize(); ++i){
      phi = Phi[i];
      break;
    }
    h_Phi->Fill(phi);
    //myfile << phi << " ";

    myfile << " \n" ;
  //  cout << endl;
  }

  myfile.close();

  TFile *f_out = new TFile(Form("out/hist_%s",inputfile.c_str()),"RECREATE");
  h_chN->Write();
  h_ch_dR->Write();
  h_ch_dEta->Write();
  h_ch_dPhi->Write();
  h_chIso->Write();
  h_nhN->Write();
  h_nh_dR->Write();
  h_nh_dEta->Write();
  h_nh_dPhi->Write();
  h_nhIso->Write();
  h_phN->Write();
  h_ph_dR->Write();
  h_ph_dEta->Write();
  h_ph_dPhi->Write();
  h_phIso->Write();
  h_pT->Write();
  h_Eta->Write();
  h_Phi->Write();

}
