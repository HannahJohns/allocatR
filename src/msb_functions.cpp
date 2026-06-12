#include <Rcpp.h>
// #include <R.h>

using namespace Rcpp;


struct Intervention_msb{

  int biasDirection = -1;
  std::vector<int> vote;
  std::vector<double> testStatistic;

  Intervention_msb(int nCov){
    vote.resize(nCov,-1);
    testStatistic.resize(nCov,-1);
  }

};

// There is probably a better way to handle this, but it corrects faults between
// R and CPP implementations
#define MSB_EPSILON 0.000000001

// We can't repeat the trick we used for csmsb where we store calculations between runs.
// We can get the sum of values in each group (to get mean/counts) but can't track
// sum of deviations from mean (to get SD) because the mean updates every time.
// Inputs are therefore more simple here than they are for csmsb because optimisation options
// will only work for edge cases

Intervention_msb msb_interventionDecision(
    std::vector<int> &dataType,
    std::vector<int> &nClasses,
    std::vector<std::vector<double>> &covariates,
    std::vector<int> &arm,
    int n0,
    int n1,
    std::vector<double> threshold,
    double parms
){


//   Rprintf("\n\nStart intervention:\n");
//
//   Rprintf("Data types are: ");
//   for(int i=0; i< dataType.size(); i++){
//     Rprintf("%d ",dataType[i]);
//   }
//
//   Rprintf("\nNumber of classes are: ");
//   for(int i=0; i< nClasses.size(); i++){
//     Rprintf("%d ",nClasses[i]);
//   }
//
//   Rprintf("\nCovariates are: \n");
//   for(int i=0;  i< covariates.size(); i++){
//     for(int j = 0; j < n0+n1; j++){
//       Rprintf("%f ", covariates[i][j]);
//     }
//     Rprintf("\n");
//   }
//
//   Rprintf("\nArms are: \n");
//   for(int i=0; i< n0+n1; i++){
//     Rprintf("%d ",arm[i]);
//   }
//
//
//   Rprintf("\nThresholds are: \n");
//   for(int i=0; i< threshold.size(); i++){
//     Rprintf("%f ",threshold[i]);
//   }
//
//   Rprintf("\nn0 is %d n1 is %d parms is %f\n", n0, n1, parms);

  int nCovs = dataType.size();

  Intervention_msb out(nCovs);

  for(int cov_i=0; cov_i<nCovs; cov_i++)
  {
    // Behaviour here depends on if we're using nominal or continuous variables
    if(dataType[cov_i] == 2)
    {

      // // Nominal - use Chisq test
      //
      // Initialise xtabs and other values we need
      std::vector<std::vector<int>> O;
      std::vector<int> E_arm(2,0);
      std::vector<int> E_cov(nClasses[cov_i],0);

      for(int arm_i=0;arm_i<2;arm_i++)
      {
        std::vector<int> tmp_row(nClasses[cov_i],0);
        O.push_back(tmp_row);
      }

      // Get crosstab of observed values and count row/column totals
      for(int i=0;i<n0+n1;i++)
      {
        int this_cov_index = std::floor(covariates[cov_i][i]+0.5);
        O[arm[i]][this_cov_index]++;

        E_arm[arm[i]]++;
        E_cov[this_cov_index]++;
      }

      // Calculate chisq test statistic

      // Rprintf("\nObserved counts are:\n");
      // for(int arm_i=0;arm_i<2;arm_i++)
      // {
      //   for(int j=0;j<nClasses[cov_i];j++)
      //   {
      //     Rprintf("%d ",O[arm_i][j]);
      //   }
      //   Rprintf("\n");
      // }


      // Rprintf("\nExpected counts are:\n");
      // for(int arm_i=0;arm_i<2;arm_i++)
      // {
      //   for(int j=0;j<nClasses[cov_i];j++)
      //   {
      //     double E = E_arm[arm_i]*E_cov[j]/((double) n0+n1);
      //
      //     Rprintf("%f ",E);
      //   }
      //   Rprintf("\n");
      // }


      out.testStatistic[cov_i] = 0;
      bool invalidTestStat = false;
      for(int arm_i=0;arm_i<2;arm_i++)
      {
        for(int j=0;j<nClasses[cov_i];j++)
        {
          double E = E_arm[arm_i]*E_cov[j]/((double) n0+n1);
          if(E < MSB_EPSILON){
            invalidTestStat = true;
            break;
          }
          out.testStatistic[cov_i] += (std::pow(O[arm_i][j] - E, 2) )/( E );
        }

        if(invalidTestStat) break;

      }

      if(invalidTestStat) out.testStatistic[cov_i] = -9999;

      // Calculate vote. Assume that all votes are counted here.
      // We'll do a last pass of this at the final step when tallying.

      int currentCov = std::floor(covariates[cov_i][n0+n1]+0.5);

      // Rprintf("\ncurrentCov is %d\n", currentCov);

      // Rprintf("\nO[0] is %f E[0] is %f, O[1] is %f E[1] is %f\n",
      //         (double) O[0][currentCov],
      //         E_arm[0]*E_cov[currentCov]/((double) n0+n1),
      //         (double) O[1][currentCov],
      //         E_arm[1]*E_cov[currentCov]/((double) n0+n1)
      //         );

      if((double) O[0][currentCov] + MSB_EPSILON < E_arm[0]*E_cov[currentCov]/((double) n0+n1)){
        out.vote[cov_i] = 0;
      } else if((double) O[1][currentCov] + MSB_EPSILON < E_arm[1]*E_cov[currentCov]/((double) n0+n1)){
        out.vote[cov_i] = 1;
      }

      // Rprintf("Vote = %d.\n", out.vote[cov_i]);

    } else if (dataType[cov_i] == 1){

      // Ordinal - not supported

      stop("Ordinal variables are not supported by MSB");
      return out;

    } else if (dataType[cov_i] == 0){



      // Continuous - use t-test

      double m0, m1, s0, s1;
      m0 = m1 = s0 = s1 = 0;

      // Rprintf("\n\nCovariate %d\n", cov_i);
      // Calculate means for both arms
      for(int i=0;i<n0+n1;i++)
      {

        // Rprintf("%d %f\n",arm[i],covariates[cov_i][i]);

        if(arm[i]==0){
          m0 += covariates[cov_i][i];
        } else {
          m1 += covariates[cov_i][i];
        }
      }

      // Rprintf("sum0 = %f, sum1=%f. n0=%d, n1=%d \n", m0, m1, n0, n1);

      m0 = m0/n0;
      m1 = m1/n1;

      // Calculate standard deviations

      for(int i=0;i<n0+n1;i++)
      {
        if(arm[i]==0){
          s0 += std::pow( m0 - covariates[cov_i][i],2);
        } else {
          s1 += std::pow( m1 - covariates[cov_i][i],2);
        }
      }

      // Because s0 and s1 are variance multiplied by (n0-1) and (n1-1) respectively,
      // we can go directly do calculating pooled standard error here
      double sp = (s0 + s1) / ((double) n0+n1-2);
      sp = std::sqrt(sp);

      // The above shortcut leads to errors if n0 or n1 is less than 2 as we're dividing by zero.
      // If this happens, replace testStatistic with a dummy value
      if(n0 < 2 || n1 < 2){
        out.testStatistic[cov_i] = -9999; // Test statistic thresholds ask if we're above threshold so this will work
      } else {
        out.testStatistic[cov_i] = std::abs(m0-m1)/(sp * std::sqrt( 1/((double) n0) + 1/((double) n1) ));
      }



      // Calculate vote. Assume that all votes are counted here.
      // We'll do a last pass of this at the final step when tallying.

      // Rprintf("m0 = %f, m1=%f, s0=%f, s1 = %f. New record %f.\n", m0, m1, s0, s1, covariates[cov_i][n0+n1]);

      if( (m0 < m1 && covariates[cov_i][n0+n1] - MSB_EPSILON > m1) ||
          (m1 < m0 && covariates[cov_i][n0+n1] + MSB_EPSILON < m1) ){

            // Vote 0 in either of these situations:
            //         m0----m1----(*)
            //   (*)---m1----m0
            out.vote[cov_i] = 0;

      } else if( (m0 < m1 && covariates[cov_i][n0+n1] + MSB_EPSILON < m0) ||
                 (m1 < m0 && covariates[cov_i][n0+n1] - MSB_EPSILON > m0) ){

            // Vote 1 in either of these situations:
            //   (*)----m0----m1
            //          m1----m0----(*)
            out.vote[cov_i] = 1;
      }

      // Rprintf("Vote = %d.\n", out.vote[cov_i]);

    } else {
      stop("Data type not recognised. Must be either numeric or factor");
      return out;

    }
  }

  //Imbalances are calculated, make intervention decision.
  // MSB works on a simple voting procedure

  int tally = 0;
  for(int cov_i=0; cov_i<nCovs; cov_i++)
  {
    if(out.testStatistic[cov_i] >= threshold[cov_i]){

      if(out.vote[cov_i]==0){
        tally--;
      } else if(out.vote[cov_i]==1){
        tally++;
      }
    }
  }

  if(tally>0)
  {
    out.biasDirection=1;
  }
  else if(tally<0)
  {
    out.biasDirection=0;
  }

  // Rprintf("BD is %d\n",out.biasDirection);

  return out;

}



// [[Rcpp::export]]
List msb_get_decision(int nObs,
                  IntegerVector dataType_in,
                  IntegerVector nClasses_in,
                  NumericMatrix data,
                  IntegerVector arm_in,
                  NumericVector threshold_in,
                  double parms
)
{

  if(data.nrow() <= nObs) Rprintf("PROVIDE THE NEW PATIENT AS WELL YA DINGUS\n");
  // Take a dataset and vote on the arm that a new patient should be assigned to
  // Should return the same values as doing this at the R level would return

  // Same structure as will be in the new setup.
  // All covariates are pre-generated in R and fed to C++

  // In reality we know where we're up to because we'll be on row n0+n1
  // (zero indexed)

  // For the purposes of this, however, we'll use nObs to indicate which row
  // we need to go to

  int nCov = dataType_in.length();

  if(threshold_in.size()<nCov) Rprintf("Provide thresholds for all variables!\n");

  std::vector<int> dataType(nCov);
  std::vector<int> nClasses(nCov);
  std::vector<double> threshold(nCov);
  std::vector<std::vector<double>> covariates;
  std::vector<int> arm(nObs,0);
  int n0=0;
  int n1=0;

  // Initialise vectors
  for(int cov_i=0;cov_i<nCov;cov_i++)
  {
    dataType[cov_i] = dataType_in[cov_i];
    nClasses[cov_i] = nClasses_in[cov_i];
    threshold[cov_i] = threshold_in[cov_i];

    std::vector<double> thisColumn(nObs+1,-9999);
    covariates.push_back(thisColumn);
  }

  // Take existing data from R and feed it into our std::vector nonsense format.
  for(int i=0;i<nObs;i++)
  {
    arm[i] = arm_in(i);
    if(arm[i]==0)
    {
      n0++;
    }
    else
    {
      n1++;
    }

    for(int cov_i=0;cov_i<nCov;cov_i++)
    {
      covariates[cov_i][i] = data(i,cov_i);
    }
  } // end for i

  // Do the same for the yet-to-be randomised patient
  for(int cov_i=0;cov_i<nCov;cov_i++)
  {
    covariates[cov_i][n0+n1] = data(n0+n1,cov_i);
  }

  // Take all of this nonsense and get a decision
  Intervention_msb intervention = msb_interventionDecision(dataType,
                                                   nClasses,
                                                   covariates,
                                                   arm,
                                                   n0,
                                                   n1,
                                                   threshold,
                                                   parms
  );


  double probability =-1;

  switch(intervention.biasDirection){
    case -1:
      probability = 0.5;
      break;
    case 0:
      probability = parms;
      break;
    case 1:
      probability = 1-parms;
      break;
  }

  // Rprintf("Probability = %f\n", probability);

  List out = List::create(
    Named("biasDirection") = intervention.biasDirection,
    Named("probability") = probability,
    Named("vote") = wrap(intervention.vote),
    Named("testStatistic") = wrap(intervention.testStatistic)
  );

  return out;
}




// [[Rcpp::export]]
List simulate_msb(IntegerVector dataType_in,
                  IntegerVector nClasses_in,
                  NumericMatrix data,
                  NumericVector roll_in,
                  NumericMatrix threshold_in,
                  double parms
){

  int nObs = data.nrow();

  int nCov = dataType_in.length();

  std::vector<int> dataType(nCov);
  std::vector<int> nClasses(nCov);
  std::vector<double> threshold(nCov);
  std::vector<double> roll(nObs);

  std::vector<std::vector<double>> covariates;
  std::vector<int> arm(nObs,0);

  for(int i=0;i<nObs;i++)
  {
    roll[i] = roll_in(i);
  }

  // Initialise vectors
  for(int cov_i=0;cov_i<nCov;cov_i++)
  {
    dataType[cov_i] = dataType_in(cov_i);
    nClasses[cov_i] = nClasses_in(cov_i);

    std::vector<double> thisColumn(nObs,-9999);
    covariates.push_back(thisColumn);
  }

  // Randomisation outputs are defined here

  // Vote placed by each covariate at each step
  IntegerMatrix vote(nObs,nCov);
  NumericMatrix testStatistic(nObs,nCov);



  // Randomisation probability (in favour of treatment)
  NumericVector randProb(nObs);

  // Load in the first patient
  for(int cov_i=0;cov_i<nCov;cov_i++)
  {
    covariates[cov_i][0] = data(0,cov_i);
  }

  int n0=0;
  int n1=0;


  // Main loop here
  for(int i=0;i<nObs;i++)
  {

    // Rprintf("\n--------------------\nRecord %d\n----------------\n",i);

    // Rprintf("\n\nRecord %d begin \n\n", i);

    // Get this set of thresholds
    for(int cov_i=0;cov_i<nCov;cov_i++)
    {
      threshold[cov_i] = threshold_in(i,cov_i);
    }

    // Randomise the new patient
    Intervention_msb intervention = msb_interventionDecision(dataType,
                                                     nClasses,
                                                     covariates,
                                                     arm,
                                                     n0,
                                                     n1,
                                                     threshold,
                                                     parms
    );

    // Rprintf("\ndecision reached, bias direction is %d\n",intervention.biasDirection);

    // Log randomisation and update current file

    for(int cov_i=0;cov_i<nCov;cov_i++)
    {
      vote(i,cov_i) = intervention.vote[cov_i];
      testStatistic(i,cov_i) = intervention.testStatistic[cov_i];
    }

    randProb(i) =-1;

    switch(intervention.biasDirection){
    case -1:
      randProb(i) = 0.5;
      break;
    case 0:
      randProb(i) = parms;
      break;
    case 1:
      randProb(i) = 1-parms;
      break;
    }

    // Rprintf("Probability = %f\n", randProb(i));

    int newArm;
    if(roll[i] < randProb(i))
    {
      newArm=0;
    }
    else
    {
      newArm=1;
    }

    // Rprintf("\n\n Rand prob = %f, roll = %f, arm = %d \n", randProb(i), roll[i], newArm);

    arm[i] = newArm;
    if(arm[i]==0)
    {
      n0++;
    }
    else
    {
      n1++;
    }

    // Add the next patient to the data frame should this not be the final randomisation
    if(i<(nObs-1))
    {
      for(int cov_i=0;cov_i<nCov;cov_i++)
      {
        covariates[cov_i][i+1] = data(i+1,cov_i);
      }
    }

  } // end for i

  List out = List::create(
    Named("arm")=arm,
    Named("vote")=vote,
    Named("testStatistic")=testStatistic,
    Named("randProb")=randProb
  );

  return out;

}

