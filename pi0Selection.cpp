//==============
// pi0 selection
//==============
for(std::vector<Mdst_pi0>::iterator it_pi0all=mdstpi0mgr.begin(); it_pi0all!=mdstpi0mgr.end(); ++it_pi0all )
{
// Number of all Mdst_pi0.
  H_nofpi0[0]->accumulate( (float)mdstpi0mgr.count(), 1.0 );
  Mdst_pi0& pi0all = *it_pi0all;

// Check each photon's momentum vector.
  Hep3Vector gamma0_3v( pi0all.gamma(0).px(), pi0all.gamma(0).py(), pi0all.gamma(0).pz() );
  Hep3Vector gamma1_3v( pi0all.gamma(1).px(), pi0all.gamma(1).py(), pi0all.gamma(1).pz() );

// Apply energy threshold.
  if( gamma_tight(gamma0_3v)==0 && gamma_tight(gamma1_3v)==0 )
  {
    Particle pi0_tight( pi0all );
    vect_pi0.push_back( pi0_tight );
  }
}

