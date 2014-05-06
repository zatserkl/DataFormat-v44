{
// gROOT->Macro(gSystem->ExpandPathName("$(HOME)/macros/rootlogon.C"));
// cout<< "*-- Local rootlogon" << endl;

cout<< "Load DataFormat.C+" <<endl;
gROOT->LoadMacro("DataFormat.C+");

cout<< "Load processRun.C+" <<endl;
gROOT->LoadMacro("processRun.C+");
}
