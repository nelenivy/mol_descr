#!/usr/bin/perl -w
#
use strict;
use MdmDiscoveryScript;

my $document = DiscoveryScript::Open(
    {
        Path      => 'E:\_disser\Diplom\Sets\Activity\glik110.sdf'
        
    }
);

my $arrayOfMolecules = $document->Molecules;

for (my $ind = 0; $ind < $arrayOfMolecules->Count(); $ind++)
{
	$arrayOfMolecules->Item($ind)->SetInvisible();
}
# get all the atoms in a molecule
for (my $ind = 0; $ind < $arrayOfMolecules->Count(); $ind++)
{
	my $molecule = $arrayOfMolecules->Item($ind);
	$molecule->SetVisible();
	my $arrayOfAtoms = $molecule->Atoms; 
	my $surface =
	  $document->CreateSolidSurface( $arrayOfAtoms, Mdm::surfaceStyleSolvent, False,
	    Mdm::surfaceColorByElectrostaticPotential, 1.4 );	
	$document->Save("E:\\srf$ind.wrl" ,'wrl' ); 
	$surface->SetInvisible();
	$molecule->SetInvisible();
	$surface = undef;
	$molecule = undef;
	$arrayOfAtoms = undef;
}
     