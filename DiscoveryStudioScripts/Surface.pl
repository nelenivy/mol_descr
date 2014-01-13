#!/usr/bin/perl -w
#
use strict;
use MdmDiscoveryScript;
my $i=1;
my $a='dgh';
 for ($i = 430; $i <= 430; $i++) {

my $document = DiscoveryScript::Open(
    {
        Path      => join('','D:\Matlab\3d\Known_',$i,'.mol2'),
        FormatType => "mol2"
    }
);
my $allAtoms = $document->Atoms;
my $surface =
  $document->CreateSolidSurface( $allAtoms, Mdm::surfaceStyleSolvent, False,
    Mdm::surfaceColorByElectrostaticPotential, 1.4 );
$document->Save(join('','D:\Matlab\3d\Known_',$i,'.wrl') ,'wrl' ); 
$document->Close();
}