#include "TCanvas.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TPad.h"

#include "Interpolenta/Interpolenta"
#include "Interpolenta/Utils"

#include <array>
#include <cmath>
#include <iostream>
#include <vector>

const double MY_PI = std::asin(1.0) * 2;

int main() {

  TCanvas c1("c1", "");

  {
    TGraph sing_g;

    TGraph sing_g_rare;

    int npoints = 1000;
    double periods = 2;

    std::vector<std::array<double, 2>> RBFs;

    int j = 0;
    for (int i = 0; i < npoints; ++i) {
      double x = (double(i) * (MY_PI * 2.0 * periods / double(npoints))) -
                 (MY_PI * periods);
      sing_g.SetPoint(i, x, std::sin(x));

      if (i && !(i % 50)) {
        sing_g_rare.SetPoint(j++, x, std::sin(x));
        RBFs.push_back({x, std::sin(x)});
      }
    }

    Eigen::MatrixXd KnotMatrix = to_knots(RBFs);

    sing_g_rare.SetMarkerColor(kRed);
    sing_g_rare.SetMarkerStyle(20);
    sing_g_rare.SetMarkerSize(1);

    sing_g.Draw("APL");
    sing_g_rare.Draw("P");

    double shape = 2;

    Interpolenta interpol(KnotMatrix);

    TGraph Interpolated_g;
    Eigen::Matrix<double, 1, 1> X;

    for (int i = 0; i < npoints; ++i) {

      X(0) = (double(i) * (MY_PI * 2.0 * periods / double(npoints))) -
             (MY_PI * periods);

      Interpolated_g.SetPoint(i, X(0), interpol.Eval(X));
    }

    Interpolated_g.SetLineColor(kBlue);
    Interpolated_g.Draw("L");

    c1.Print("canv1d.pdf");
  }

  {
    c1.Clear();
    c1.cd();
    TPad p1("p1", "", 0, 0, 0.5, 0.5);
    p1.AppendPad();

    c1.cd();
    TPad p2("p2", "", 0.5, 0.5, 1, 1);
    p2.AppendPad();

    c1.cd();
    TPad p3("p3", "", 0, 0.5, 0.5, 1);
    p3.AppendPad();

    c1.cd();
    TPad p4("p4", "", 0.5, 0, 1, 0.5);
    p4.AppendPad();

    TGraph2D sing_g;

    TGraph sing_g_rare;
    TGraph2D sing2_g_rare;

    int npoints = 100;
    double periods = 2;

    std::vector<std::array<double, 3>> RBFs;

    size_t n_data_points = 0;
    size_t n_points = 0;
    for (int i = 0; i < npoints; ++i) {
      for (int j = 0; j < npoints; ++j) {
        double x = (double(i) * (MY_PI * 2.0 * periods / double(npoints))) -
                   (MY_PI * periods);
        double y = (double(j) * (MY_PI * 2.0 * periods / double(npoints))) -
                   (MY_PI * periods);
        double val = std::sin(x) + std::sin(y);
        sing_g.SetPoint(n_points++, x, y, val);

        if (i && !(i % 10) && j && !(j % 10)) {
          sing_g_rare.SetPoint(n_data_points, x, y);
          sing2_g_rare.SetPoint(n_data_points++, x, y, val);
          RBFs.push_back({x, y, val});
        }
      }
    }

    sing_g_rare.SetMarkerColor(kRed);
    sing_g_rare.SetMarkerStyle(20);
    sing_g_rare.SetMarkerSize(1);

    sing2_g_rare.SetMarkerColor(kRed);
    sing2_g_rare.SetMarkerStyle(20);
    sing2_g_rare.SetMarkerSize(1);

    p1.cd();
    sing_g.Draw("COLZ");
    sing_g_rare.Draw("P");

    p2.cd();
    sing2_g_rare.Draw("P0");

    Eigen::MatrixXd KnotMatrix = to_knots(RBFs);
    Interpolenta interpol(KnotMatrix);
    Eigen::Matrix<double, 2, 1> X;

    TGraph2D Interpolated_g;

    int npoints_check = 0;
    for (int i = 0; i < npoints; ++i) {
      for (int j = 0; j < npoints; ++j) {
        X(0) = (double(i) * (MY_PI * 2.0 * periods / double(npoints))) -
               (MY_PI * periods);
        X(1) = (double(j) * (MY_PI * 2.0 * periods / double(npoints))) -
               (MY_PI * periods);

        Interpolated_g.SetPoint(npoints_check++, X(0), X(1), interpol.Eval(X));
      }
    }

    Interpolated_g.SetMarkerStyle(20);
    Interpolated_g.SetMarkerSize(1);

    p3.cd();
    Interpolated_g.Draw("PCOL");

    p4.cd();
    Interpolated_g.Draw("COLZ");

    c1.Print("canv2d.pdf");
  }

  return 1;
}