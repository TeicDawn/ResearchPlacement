#include "SetNiceStyle.C"

void MapChi2Pois(string PID = "track") { //Segun la label, track, shower o ninguna de ellas, os hace esto para tracks, showers y middles

  TFile * f1 = new TFile("Events_1y_NO_Std_upto20.root"); //Aquí va el nombre de vuestro archivo, si no está en este path, dais el path entero
  auto hStdt = (TH2D * ) f1 -> Get("h_tracks_all_yx") -> Clone(""); //En los archivos que os mandaremos se llamará así el histograma en cos-theta/E
  auto hStds = (TH2D * ) f1 -> Get("h_showers_all_yx") -> Clone(""); //Data
  auto hStdm = (TH2D * ) f1 -> Get("h_middles_all_yx") -> Clone(""); //Data
  TFile * f2 = new TFile("Events_1y_NO_Decay1e-5_upto20.root"); //El de arriba será para Data, el de abajo para Model
  auto hNSIt = (TH2D * ) f2 -> Get("h_tracks_all_yx") -> Clone(""); //Este para Model
  auto hNSIs = (TH2D * ) f2 -> Get("h_showers_all_yx") -> Clone(""); //Model
  auto hNSIm = (TH2D * ) f2 -> Get("h_middles_all_yx") -> Clone(""); //Model
  auto c = new TCanvas();

  hStdt -> Print();
  hNSIt -> Print();

  //Para el pseudo-experimento, los histogramas se cogen así, son en 3D (Bjorken-y)
  auto hToy = (TH3D * ) f2 -> Get("h_tracks_toy_all") -> Clone("");
  //Pero solo tienen un bin en Bj, así que simplemente GetBinContent(i,j,1) y todo igual.

  //Esto para quedarnos solo con upgoing events

  hNSIt -> GetYaxis() -> SetRangeUser(-1, 0);
  hNSIs -> GetYaxis() -> SetRangeUser(-1, 0);
  hNSIm -> GetYaxis() -> SetRangeUser(-1, 0);
  hStdt -> GetYaxis() -> SetRangeUser(-1, 0);
  hStds -> GetYaxis() -> SetRangeUser(-1, 0);
  hStdm -> GetYaxis() -> SetRangeUser(-1, 0);

  //Tomamos los n bins en x e y
  int nbinsx = hStdt -> GetNbinsX();
  int nbinsy = hStdt -> GetNbinsY();
  hStdt -> Print();
  hNSIt -> Print();
  double sens = 0; //Aquí computaremos la chi2
  if (PID == "track") {
    for (int i = 1; i <= nbinsx; i++) {
      for (int j = 1; j <= nbinsy; j++) {

        double std = hStdt -> GetBinContent(i, j); //Este será para Data, ignorad los nombres de std, str, hStd, hNSI
        double str = hNSIt -> GetBinContent(i, j); //Este para Model
        cout << std << " " << str << endl;

        double sigma = 0;

        if (std > 0) sigma = 2*(str - std + str*log(str/std)); //Esta es la poissoniana
        if (std == 0) sigma = 0;
        if (sigma > 100 or sigma < -100) sigma = 0;
        //Cambiamos cada bin del histograma por la chi2
        hNSIt -> SetBinContent(i, j, sigma * fabs(sigma)); //Hacemos sigma*abs(sigma) para ver si es exceso/defecto de eventos

        sens += sigma; //Elevamos al cuadrado la sigma y sumamos para tener la total chi2
        cout << sigma * fabs(sigma) << endl;
      }
    }

    cout << "Total Chi2 = " << sens << endl;

    SetHist(hNSIt);

    hNSIt -> SetTitle(";Energy [GeV]; cos#theta_{z}; S_{ij} *|S_{ij}|");

    SetNiceStyle();
    SetNiceTempPalette();

    double max = hNSIt -> GetMaximum();
    double min = hNSIt -> GetMinimum();

    if (max < fabs(min)) max = fabs(min);
    cout << "Maximum chi2 in a bin is: " << max << endl;
    //Esto es para que se vea la chi2 con colores y 0 sea blanco, rojo exceso, azul defecto
    hNSIt -> GetZaxis() -> SetRangeUser(-max, max);
    hNSIt -> GetZaxis() -> SetNdivisions(6);
    hNSIt -> GetYaxis() -> SetRangeUser(-1, 0);
    hNSIt -> GetYaxis() -> SetNdivisions(506); //TStyle::SetNdivisons(506, "xyz") not working
    hNSIt -> Draw("colz"); //Dibujar mapa de colores

    hNSIt -> GetXaxis() -> SetMoreLogLabels();

    gPad -> SetLogx();

    gPad -> SetTopMargin(0.09);
    //Para escribir en mitad
    MiscText(0.5, 0.95, 0.05, "Tracks (Stat. only #chi^{2}) ", kBlack, true);

    //Lo guarda como png pero dependiendo de cómo lo llaméis a veces la calidad no será la deseada, podéis comentar la línea y guardarlo como vosotros queráis.
    //c->SaveAs("MapChi2tracks.png");

  } else if (PID == "shower") { //Ya es repetir para showers, podría haber hecho el código más compacto pero... requiere dolores de cabeza en c++

    for (int i = 1; i <= nbinsx; i++) {
      for (int j = 1; j <= nbinsy; j++) {

        double std = hStds -> GetBinContent(i, j);
        double str = hNSIs -> GetBinContent(i, j);

        double sigma = 0;

        if (std > 0) sigma = (str - std) / sqrt(std);

        hNSIs -> SetBinContent(i, j, sigma * fabs(sigma));

        sens += pow(sigma, 2);

      }
    }

    cout << "Total Chi2 = " << sens << endl;

    SetHist(hNSIs);

    hNSIs -> SetTitle(";Energy [GeV]; cos#theta_{z}; S_{ij} *|S_{ij}|");

    SetNiceStyle();
    SetNiceTempPalette();

    double max = hNSIs -> GetMaximum();
    double min = hNSIs -> GetMinimum();

    if (max < fabs(min)) max = fabs(min);
    cout << "Maximum chi2 in a bin is: " << max << endl;
    //max = 0.02;

    cout << max << endl;
    hNSIs -> GetZaxis() -> SetRangeUser(-max, max);
    hNSIs -> GetZaxis() -> SetNdivisions(6);
    hNSIs -> GetYaxis() -> SetRangeUser(-1, 0);
    hNSIs -> GetYaxis() -> SetNdivisions(506); //TStyle::SetNdivisons(506, "xyz") not working
    hNSIs -> Draw("colz");

    hNSIs -> GetXaxis() -> SetMoreLogLabels();

    gPad -> SetLogx();

    gPad -> SetTopMargin(0.09);

    MiscText(0.5, 0.95, 0.05, "Showers (Stat. only #chi^{2}) ", kBlack, true);

    //c->SaveAs("MapChi2showers.png");

  } else {
    for (int i = 1; i <= nbinsx; i++) {
      for (int j = 1; j <= nbinsy; j++) {

        double std = hStdm -> GetBinContent(i, j);
        double str = hNSIm -> GetBinContent(i, j);

        double sigma = 0;

        if (std > 0) sigma = (str - std) / sqrt(std);

        hNSIm -> SetBinContent(i, j, sigma * fabs(sigma));

        sens += pow(sigma, 2);

      }
    }

    cout << "Total Chi2 = " << sens << endl;

    SetHist(hNSIm);

    hNSIm -> SetTitle(";Energy [GeV]; cos#theta_{z}; S_{ij} *|S_{ij}|");

    SetNiceStyle();
    SetNiceTempPalette();

    double max = hNSIm -> GetMaximum();
    double min = hNSIm -> GetMinimum();

    if (max < fabs(min)) max = fabs(min);
    cout << "Maximum chi2 in a bin is: " << max << endl;
    //max = 0.02;

    hNSIm -> GetZaxis() -> SetRangeUser(-max, max);
    hNSIm -> GetZaxis() -> SetNdivisions(6);
    hNSIm -> GetYaxis() -> SetRangeUser(-1, 0);
    hNSIm -> GetYaxis() -> SetNdivisions(506); //TStyle::SetNdivisons(506, "xyz") not working
    hNSIm -> Draw("colz");

    hNSIm -> GetXaxis() -> SetMoreLogLabels();

    gPad -> SetLogx();

    gPad -> SetTopMargin(0.09);

    MiscText(0.5, 0.95, 0.05, "Middles (Stat. only #chi^{2}) ", kBlack, true);

    //c->SaveAs("MapChi2middles.png");

  }
}
