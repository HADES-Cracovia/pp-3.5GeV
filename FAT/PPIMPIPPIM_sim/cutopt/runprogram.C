void runprogram()
{
  gROOT->ProcessLine(".L ppimpippim.C++");
  ppimpippim l;
  l.Loop();
}
