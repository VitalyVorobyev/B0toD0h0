//
//  mystat.h
//
//  $Id$
//

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

#define MYSTAT_NCOUNT 10

class mystat {
 public:
  mystat(const char *name = "mystat");
  ~mystat() {}
  void begin_run(const char *str=0);
  void event();
  void count(int i=0, int step=1);
  void term();
  void everyevent(int every=1) { m_every = every; }
  void countname(int i, const char *name);
  
 private:
  void printline(time_t t0, int nevent, const char *add=0);
  char *m_name;
  char *m_countname[MYSTAT_NCOUNT];
  int *m_shmp;
  int m_every; // default=every 1000 event
};

#if defined(BELLE_NAMESPACE)
}
#endif
