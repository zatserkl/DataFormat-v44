{
gROOT->Macro(gSystem->ExpandPathName("$(HOME)/macros/rootlogon.C"));
cout<< "*-- Local rootlogon" << endl;

if (gClassTable->GetID("PCTEvent") < 0) {
    cout<< "Load DataFormat.C+" <<endl;
    gROOT->LoadMacro("DataFormat.C+");
}

if (gClassTable->GetID("RecoEvent") <= 0) {
    cout<< "Load recoRun.C+" <<endl;
    gROOT->LoadMacro("recoRun.C+");
}

//cout<< "Load online.C+" <<endl;
//gROOT->LoadMacro("online.C+");
}
