#!/usr/bin/perl -w
#
use strict;
use MdmDiscoveryScript;
my $i=1;
my $a='dgh';
 for ($i = 1; $i <= 449; $i++) {

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
my $molecule = $document->Molecules->Item(0);
$document->CalculateCharges();
my $j=open(my $fh, '>', join('','D:\Matlab\3d\Known_',$i,'.txt')) ;
print $fh  join('',$molecule->AtomCount,"\n");
foreach my $atom (@$allAtoms)
{
    print $fh join(' ',$atom->XYZ->X,$atom->XYZ->Y,$atom->XYZ->Z,$atom->Charge,"\n");
}
close $fh;
$document->Close();
}