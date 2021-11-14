#include "EnergyBackground_Level.hh"


EnergyBackground_Level::EnergyBackground_Level()
{
  Exb=0*keV;
  Eb=0*keV;
  fracb=0;
}

EnergyBackground_Level::~EnergyBackground_Level()
{
  ;
}
/*--------------------------------------------------------*/
void EnergyBackground_Level::ReportBackground()
{
  G4cout<<"----> Excitation energy       "<<Exb/keV<<" [keV]"<<G4endl;
  G4cout<<"   => Decay energy            "<<Eb/keV<<" [keV]"<<G4endl;
  G4cout<<"   => Direct feeding fraction "<<fracb<<" [%]"<<G4endl;
}

