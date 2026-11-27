#include <library/preparation/Preparation.h>
#include <library/preparation/AssayPrediction.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <iostream>
int main()
{
  int rejected = 0;
  try { OpenMS::OpenDIA::Preparation::AssayPredictor predictor; }
  catch (const OpenMS::Exception::Precondition&) { ++rejected; }
  try { (void)OpenMS::OpenDIA::Preparation::predictDecoys({}); }
  catch (const OpenMS::Exception::Precondition&) { ++rejected; }
  if (rejected != 2) return 1;
  std::cout << "PASS: ONNX-disabled predictor and preparation fail with Precondition\n";
}
