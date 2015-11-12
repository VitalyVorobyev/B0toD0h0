
	      //////////////////////////////////////////
	      //                                      //
	      //  rooksfw.h                           //
	      //                                      //
	      //  ROOT version of the k_sfw           //
	      //  improved Super-Fox-Wolfram moments  //
	      //  to be used with ksfwmoments.{cc,h}  //
	      //                                      //
	      //  M. Nakao                            //
	      //                                      //
	      //////////////////////////////////////////

// Versions
// 2006111301  first version
// 2007062700  float version
// 2007083000  fix (float double?)

#include <TH1.h>
#include <TH2.h>

// ----------------------------------------------------------------------
// rooksfw
// ----------------------------------------------------------------------
class rooksfw {
public:
  rooksfw() { clear(); }
  rooksfw(const rooksfw &X);
  ~rooksfw() { if (m_name) delete m_name; }
  rooksfw(const char *name, const double alpha[7][17],
	  const double sigpdf[7][8], const double bkgpdf[7][8]);

  void clear();
  void input(const double mm2, const double var[17]);
  void input(const double mm2, const float var[17]);
  void input(const double mm2, const double et,
	     const float hso[3][5], const float hoo[]);
  void input(const double mm2, const double et,
	     const double hso[3][5], const double hoo[]);
  void input(const double mm2, const double et,
	     const double hso00, const double hso01, const double hso02,
	     const double hso03, const double hso04, const double hso10,
	     const double hso12, const double hso14, const double hso20,
	     const double hso22, const double hso24,
	     const double hoo0, const double hoo1, const double hoo2,
	     const double hoo3, const double hoo4);
  
  TH1F *histmean(int i) const {
    if (i < 0 || i>= 7) return 0;
    return m_histmean[i]; }
  TH2F *histcorr(int i) const {
    if (i < 0 || i>= 7) return 0;
    return m_histcorr[i]; }
  TH1F *histksfw(int i) const {
    if (i < 0 || i>= 7) return 0;
    return m_histksfw[i]; }
  TH1F *histimm2() const { return m_histimm2; }
  TH1F *hist(int i) const {
    if (i == 21) return histimm2();
    switch (i/7) {
    case 0: return histmean(i);
    case 1: return (TH1F *)histcorr(i-7);
    default: return histksfw(i-14);
    }
  }
  void Write() { for (int i=0; i<=21; i++) hist(i)->Write(); }
  
  void train();
  void fill();

  double var(int i) const { return (i < 0 || i >= 17) ? 0 : m_var[i]; }
  double mm2() const { return m_mm2; }
  int    imm2() const { return m_imm2; }
  double et() const { return m_var[0]; }
  double hso(int i, int j) const {
    switch (i) {
    case 0: return (j<0||j>=3) ? 0 : m_var[j+1];
    case 1: return (j==0) ? m_var[4] : 0;
    case 2: return (j<0||j>=3) ? 0 : m_var[j+5];
    case 3: return (j==0) ? m_var[8] : 0;
    case 4: return (j<0||j>=3) ? 0 : m_var[j+9];
    default: return 0;
    }
  }
  double hoo(int i) const { return (i<0||i>=5) ? 0 : m_var[i+12]; }
  double ksfw() const { return m_ksfw; }
  double ls() const { return m_ls; }
  double lb() const { return m_lb; }
  void mm2_correction(double c[7]) {
    for (int i=0; i<7; i++) m_mm2_correction_sig[i] = c[i];
    m_mm2_correction = 1; }
  void mm2_correction(int a) { m_mm2_correction = a; }
  int  mm2_correction() const { return m_mm2_correction; }
  
private:
  void calc(double mm2);

  char *m_name;
  int m_filled;  // to avoid double histogram entry
  int m_trained; // to avoid double training entry
  
  TH1F *m_histmean[7];
  TH2F *m_histcorr[7];
  TH1F *m_histksfw[7];
  TH1F *m_histimm2;

  double m_alpha[7][17];
  double m_sigpdf[7][8];
  double m_bkgpdf[7][8];

  int    m_imm2;
  double m_mm2;
  double m_var[17];
  
  double m_ksfw;
  double m_ls;
  double m_lb;

  static double m_mm2_correction_sig[7];
  int    m_mm2_correction;
};
