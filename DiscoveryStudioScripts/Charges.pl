#
use strict;
use MdmDiscoveryScript;
my $i=1;
my $a='dgh';
my $document = DiscoveryScript::Open(
    {
        Path      => join('','D:\Matlab\3d\Known_',$i,'.mol2'),
        FormatType => "mol2"
    }
);
my $molecule = $document->Molecules->Item(0);
my $allAtoms = $document->Atoms;
$document->CalculateCharges();
my $j=open(my $fh, '>', join('','D:\Matlab\3d\Known_',$i,'.txt')) ;
print $j;
print $fh  join('',$molecule->AtomCount,"\n");
foreach my $atom (@$allAtoms)
{
    print $fh join(' ',$atom->XYZ->X,$atom->XYZ->Y,$atom->XYZ->Z,$atom->Charge,"\n");
}
close $fh;
$document->Close();

 for ($i = 360; $i <= 449; $i++) {

$document = DiscoveryScript::Open(
    {
        Path      => join('','D:\Matlab\3d\Known_',$i,'.mol2'),
        FormatType => "mol2"
    }
);
$molecule = $document->Molecules->Item(0);
$allAtoms = $document->Atoms;
$document->CalculateCharges();
$j=open(my $fh, '>', join('','D:\Matlab\3d\Known_',$i,'.txt')) ;
print $j;
print $fh  join('',$molecule->AtomCount,"\n");
foreach my $atom (@$allAtoms)
{
    print $fh join(' ',$atom->XYZ->X,$atom->XYZ->Y,$atom->XYZ->Z,$atom->Charge,"\n");
}
close $fh;
$document->Close();
}