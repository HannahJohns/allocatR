#include <Rcpp.h>

double getProbability(double x, int method, std::vector<double> parms)
{
  double out;

  switch(method)
  {
    case 0:
      out = parms[0];
    break;
    case 1:
      out = 0.5 + parms[0] * x;
    break;
    case 2:
      out = 0.5 + parms[0] * pow(x, parms[1]); // ax^b
    break;
    case 3:
      out = 0.5 + parms[0] * exp(parms[1] * x); // ae^(bx)
    break;
    case 4:
      out = 0.5 + parms[0] * log(parms[1] * x); // aln(bx)
    break;
    case 5:
      out = parms[0]+x*parms[1];
      if(out>parms[3]) out = parms[3];
      if(out<parms[2]) out = parms[2]; // max: parms[3]; min: parms[2]; a: parms[0]; b: parms[1]
    break;
    default:
      out=0.5;
  }

  // Can't have probabilities greater than one
  if(out>1) out = 1;

  // Rprintf("%f %f\n",x, out);

  return out;
}
