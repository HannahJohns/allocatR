#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List pocock_randomisation(const IntegerMatrix& data,
                          const IntegerVector& stratum,
                          const int nStrata,
                          const NumericVector& roll,
                          const IntegerVector& nClasses,
                          double biasStrength,
                          bool balanceGroups)
{


  std::vector<int> n0(nStrata,0);
  std::vector<int> n1(nStrata,0);

  int nObs = data.nrow();
  int nVars = data.ncol();

  IntegerVector arm(nObs);
  NumericVector prob(nObs);
  IntegerVector direction(nObs);

  // Hoo boy
  // currentBalance[stratum][covariate][covariate value][arm]

  std::vector<std::vector<std::vector<std::vector<int>>>> currentBalance;

  for(int stratum_i = 0; stratum_i<nStrata; stratum_i++){

    std::vector<std::vector<std::vector<int>>> stratum_balance;

    for(int i=0;i<nVars;i++){
      std::vector<std::vector<int>> thisTable;
      for(int j=0; j<nClasses(i);j++)
      {
        std::vector<int> thisRow = {0,0};
        thisTable.push_back(thisRow);
      }
      stratum_balance.push_back(thisTable);
    }

    currentBalance.push_back(stratum_balance);
  }

  // Run randomisations
  for(int i=0; i<nObs; i++)
  {
    // Rprintf("==========================================\n");
    // Rprintf("Assigning patient:  ");
    // for(int i_cov=0;i_cov<nVars;i_cov++)
    // {
    //   Rprintf("%d ",data(i,i_cov));
    // }
    // Rprintf("\n\n");
    //
    //
    // for(int i_cov=0;i_cov<nVars;i_cov++)
    // {
    //   for(int j_class=0; j_class< nClasses(i_cov); j_class++)
    //   {
    //     for(int k_arm=0; k_arm<2; k_arm++)
    //     {
    //       Rprintf("%d ",currentBalance[i_cov][j_class][k_arm]);
    //     }
    //     if(j_class == data(i,i_cov)) Rprintf("*");
    //     Rprintf("\n");
    //   }
    //   Rprintf("\n - - - - - - \n");
    // }

    //Get current imbalance
    // negative indicates vote for arm 0
    // positive indicates vote for arm 1
    int vote_balance =0 ;

    if(balanceGroups)
    {
      if(n0[stratum(i)]>n1[stratum(i)])
      {
        vote_balance = 1;
      }
      else if(n1[stratum(i)] > n0[stratum(i)])
      {
        vote_balance = -1;
      }
    }

    for(int i_cov=0;i_cov<nVars;i_cov++)
    {
      if(currentBalance[stratum(i)][i_cov][data(i,i_cov)][0] > currentBalance[stratum(i)][i_cov][data(i,i_cov)][1])
      {
        vote_balance++;
      }
      else if (currentBalance[stratum(i)][i_cov][data(i,i_cov)][0] < currentBalance[stratum(i)][i_cov][data(i,i_cov)][1])
      {
        vote_balance--;
      }
    }
    // Rprintf("\n Current vote balance is %d\n",vote_balance);

    direction(i) = vote_balance;

    // Get control probability
    double thisProb = (vote_balance>0) ? (1-biasStrength) : ( (vote_balance<0) ? biasStrength : 0.5);

    prob(i) = thisProb;

    // Assign treatment arm
    if(roll(i) < thisProb)
    {
      arm(i) = 0;
      n0[stratum(i)]++;
    }
    else
    {
      arm(i) = 1;
      n1[stratum(i)]++;
    }

    // Rprintf("\n Control probability is %f, roll is %f\n",thisProb,roll(i));
    // Rprintf("\n Arm is %d\n",arm(i));


    // Update covariate balance
    for(int i_cov=0;i_cov<nVars;i_cov++)
    {
      currentBalance[stratum(i)][i_cov][data(i,i_cov)][arm(i)]++;
    }
  }

  List out = List::create(
    Named("arm") = arm,
    Named("prob") = prob,
    Named("direction") = direction
  );

  return out;
}
