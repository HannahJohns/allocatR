#include "biasStrengthControl.hpp"
using namespace Rcpp;

template <typename T> int sgn(T val) {
  return (T(0) < val) - (val < T(0));
}

// There is probably a better way to handle this, but it corrects faults between
// R and CPP implementations
#define CSMSB_EPSILON 0.000000001

void genodds( std::vector<double> &out,
          std::vector<int> &Nconc,
          std::vector<int> &Ndisc,
          std::vector<int> &Ntied_conc,
          std::vector<int> &Nwin,
          std::vector<int> &Nloss,
          std::vector<int> &Ntied,
          std::vector<int> &arm,
          int n0,
          int n1
          )
{
  /*
    Calculate genodds given counts of concordant/discordant/etc for every
    observation

    This is somewhat messy because we're counting concordant/discordant pairs
    not probabilities of concordance

    i.e. Rs = 1/(n0+n1) * Nconc

  */

  double Pc = 0;
  double Pd = 0;
  double Pc_null = 0;
  double Pd_null = 0;

  // Get Pc and Pd
  for(int i=0; i<n0+n1; i++)
  {
    Pc += Nconc[i] + 0.5 * Ntied_conc[i];
    Pd += Ndisc[i] + 0.5 * Ntied_conc[i];

    // We need to spread each obs across both arms

    for(int i_arm = 0 ; i_arm <2; i_arm++)
    {
      double this_Nwin_XNull =  (i_arm==0 ? 0 : n0);
      double this_Nloss_XNull = (i_arm==0 ? n1 : 0);
      double this_comparable_XNull = (i_arm==0 ? n1 : n0);

      double this_NconcNull = this_Nwin_XNull * Nwin[i] + this_Nloss_XNull * Nloss[i];
      double this_NdiscNull = this_Nwin_XNull * Nloss[i] + this_Nloss_XNull * Nwin[i];
      double this_Ntied_concNull = this_comparable_XNull * Ntied[i];

      // Scaling by N^2 returns us to probability of concordant/discordant
      // Following the maths, this pops out of the fact that we're multiplying
      // two probabilities together

      this_NconcNull = this_NconcNull * (i_arm == 0 ? n0 : n1) / ((double) (n0+n1)*(n0+n1));
      this_NdiscNull = this_NdiscNull * (i_arm == 0 ? n0 : n1) / ((double) (n0+n1)*(n0+n1));
      this_Ntied_concNull = this_Ntied_concNull * (i_arm == 0 ? n0 : n1) / ((double) (n0+n1)*(n0+n1));

      Pc_null += this_NconcNull + 0.5 * this_Ntied_concNull;
      Pd_null += this_NdiscNull + 0.5 * this_Ntied_concNull;
    }

  }

  if(Pd==0)
  {
    // Return magic number error code if genodds is undefined
    out[0]=out[1]=out[2]=-1;
    return;
  }

  Pc = Pc/( (n0+n1)*(n0+n1) );
  Pd = Pd/( (n0+n1)*(n0+n1) );
  Pc_null = Pc_null/( (double) (n0+n1)*(n0+n1));
  Pd_null = Pd_null/( (double) (n0+n1)*(n0+n1));

  // Rcout <<  Pc << "  " << Pd << std::endl;
  // Rcout <<  Pc_null << "  " << Pd_null << std::endl;

  double OR = Pc/Pd;
  double OR_null = 1;

  double SE = 0; // not actually needed and can probably be dropped
  double SE_null = 0;

  // Standard error
  for(int i=0; i<n0+n1; i++)
  {
    SE += (OR * (Ndisc[i] + 0.5 * Ntied_conc[i]) - (Nconc[i]  + 0.5 * Ntied_conc[i])) *
          (OR * (Ndisc[i] + 0.5 * Ntied_conc[i]) - (Nconc[i]  + 0.5 * Ntied_conc[i])) /
          ( (n0+n1)*(n0+n1) );


    for(int i_arm = 0 ; i_arm <2; i_arm++)
    {

      double this_Nwin_XNull =  (i_arm==0 ? 0 : n0);
      double this_Nloss_XNull = (i_arm==0 ? n1 : 0);
      double this_comparable_XNull = (i_arm==0 ? n1 : n0);

      double this_NconcNull = this_Nwin_XNull * Nwin[i] + this_Nloss_XNull * Nloss[i];
      double this_NdiscNull = this_Nwin_XNull * Nloss[i] + this_Nloss_XNull * Nwin[i];
      double this_Ntied_concNull = this_comparable_XNull * Ntied[i];

      this_NdiscNull = (this_NdiscNull + 0.5 * this_Ntied_concNull) / ((n0+n1)*(n0+n1));
      this_NconcNull = (this_NconcNull + 0.5 * this_Ntied_concNull) / ((n0+n1)*(n0+n1));

      // Extra divide by N is needed here because we're proportionally distributing every
      // patient across both arms. the Added division rescales the first part of this
      // calculation back down to where it's meant to be

      SE_null += (i_arm == 0 ? n0 : n1) *
                 (OR_null * this_NdiscNull - this_NconcNull) *
                 (OR_null * this_NdiscNull - this_NconcNull) /
                 (n0+n1);
    }
  }

  SE = 2/(Pd * (n0+n1)) * sqrt(SE);
  SE_null = 2/(Pd_null * (n0+n1)) * sqrt(SE_null);

  out[0] = OR;
  out[1] = SE;
  out[2] = SE_null;
}


// Given a bunch of inputs for genodds, make a decision about whether or not to intervene
// in randomisation

// In a perfect world where nominal variables didn't exist, we'd have a simple
// case of one vector per variable. Instead, we need to do some nonsense to get this
// to work.

// Each input is a twice-nested vector. The outer level is the covariate number
// and the next level down is the version of the covariate we care about.

// For sensible covariates, the actual vector will be accessed as Nconc[cov_i][0]
// For nominal covariates, the relevant vector will be accessed as Nconc[cov_i][cov_value]

// where cov_i is the index of the covariate and cov_value is the value of the nominal variable

// dataType -- an integer vector giving magic numbers for thge type of covariate

    // 0 for continuous (run as-is)
    // 1 for ordinal (check for rounding issues)
    // 2 nominal (dummy code it)

// nClasses -- an integer vector giving the number of cases the nominal covariate can take on


// We should also offer controls for the full family of MSB approaches

// Return an object of type Intervention, defined as follows:

struct csmsb_Intervention{
  int biasDirection;
  double imbalanceStrength;
  double biasProbability;

  double biasProbability_covariate;
  double biasProbability_group;

  std::vector<int> vote;
  std::vector<double> testStatistic;
  std::vector<std::vector<double>> imbalance; // OR, SE(OR), SEnull(OR)
  std::vector<std::vector<std::vector<double>>> imbalance_prospective;

  // Explanatory information for debug/transparency purposes

  // FLAG: GET RID OF THESE
  double BD;
  double seBD;

  csmsb_Intervention(int nCov){
    vote.resize(nCov,-1);
    testStatistic.resize(nCov,-1);
  }

};


csmsb_Intervention csmsb_interventionDecision(
                          std::vector<int> &dataType,
                          std::vector<int> &nClasses,
                          std::vector<std::vector<double>> &covariates,
                          std::vector<int> &arm,
                          int n0,
                          int n1,
                          std::vector<double> threshold,
                          std::vector<double> parms,
                          int voteMethod,
                          int biasMethod,
                          std::vector<double> weights,
                          double groupImbalanceThreshold,
                          int imbalancePriority ,
                          std::vector<std::vector<std::vector<int>>> &Nconc,
                          std::vector<std::vector<std::vector<int>>> &Ndisc,
                          std::vector<std::vector<std::vector<int>>> &Ntied_conc,
                          std::vector<std::vector<std::vector<int>>> &Nwin,
                          std::vector<std::vector<std::vector<int>>> &Nloss,
                          std::vector<std::vector<std::vector<int>>> &Ntied
                          )
{

  int nCovs = dataType.size();

  csmsb_Intervention out(nCovs);

  for(int cov_i=0; cov_i<nCovs; cov_i++)
  {

    // For nominal variables, cast to int to get the index for
    // genodds components based on the value of the most recent patient
    // Otherwise it's just the first index
    int dummyCodeIndex = (dataType[cov_i] == 2 ? std::floor(covariates[cov_i][n0+n1]+0.5) : 0);

    if(dataType[cov_i] == 2 && dummyCodeIndex >= nClasses[cov_i])
    {
      stop("Nominal variable code was greater than maximum class value");
      return out;
    }


    // Calculate prospective imbalance if assigned to either arm
    // Need to get prospective pairs in both cases, resulting in
    // yet more nested vectors

    std::vector<std::vector<int>> Nconc_prospective;
    std::vector<std::vector<int>> Ndisc_prospective;
    std::vector<std::vector<int>> Ntied_conc_prospective;
    std::vector<std::vector<int>> Nwin_prospective;
    std::vector<std::vector<int>> Nloss_prospective;
    std::vector<std::vector<int>> Ntied_prospective;
    std::vector<std::vector<int>> arm_prospective;

    for(int arm_i=0;arm_i<2;arm_i++)
    {
      std::vector<int> tmp_prospective(n0+n1+1,0);

      Nconc_prospective.push_back(tmp_prospective);
      Ndisc_prospective.push_back(tmp_prospective);
      Ntied_conc_prospective.push_back(tmp_prospective);
      Nwin_prospective.push_back(tmp_prospective);
      Nloss_prospective.push_back(tmp_prospective);
      Ntied_prospective.push_back(tmp_prospective);
      arm_prospective.push_back(tmp_prospective);
    }


    for(int arm_i=0;arm_i<2;arm_i++)
    {
      // Copy over current state, then
      // assume the new patient is assigned to group 0,
      // update concordant/discordant/etc pair counts here

      arm_prospective[arm_i][n0+n1] = arm_i;
      for(int i=0; i<(n0+n1); i++)
      {

        // Copy over the values we care about so we can update them
        Nconc_prospective[arm_i][i] = Nconc[cov_i][dummyCodeIndex][i];
        Ndisc_prospective[arm_i][i] = Ndisc[cov_i][dummyCodeIndex][i];
        Ntied_conc_prospective[arm_i][i] = Ntied_conc[cov_i][dummyCodeIndex][i];
        Nwin_prospective[arm_i][i] = Nwin[cov_i][dummyCodeIndex][i];
        Nloss_prospective[arm_i][i] = Nloss[cov_i][dummyCodeIndex][i];
        Ntied_prospective[arm_i][i] = Ntied[cov_i][dummyCodeIndex][i];
        arm_prospective[arm_i][i] = arm[i];


        // Update with the prospective assignment

        // Does the new patient win/lose on both x and y?
        int indicator_x = arm_prospective[arm_i][n0+n1] - arm_prospective[arm_i][i];

        double epsilon = 0.0001;
        int indicator_y;
        if(dataType[cov_i]==0)
        {
          indicator_y = sgn(covariates[cov_i][n0+n1] - covariates[cov_i][i]);
        }
        else if(dataType[cov_i]==1)
        {

          indicator_y = (
            (abs(covariates[cov_i][n0+n1] - covariates[cov_i][i]) < epsilon) ?
          0 :
          sgn(covariates[cov_i][n0+n1] - covariates[cov_i][i])
          );
        }
        else if(dataType[cov_i]==2)
        {
          // We know that the dummy coding for the prospecive patient
          // is 1 by definition, so all we need to do is check if it's
          // a match with each existing patient.
          // If it's a match, then we're tied. If it's not a match, then
          // the new patient wins on this dummy coding.

          indicator_y = (
            (abs(covariates[cov_i][n0+n1] - covariates[cov_i][i]) < epsilon) ? 0:1
          );

        }
        else
        {
          Rprintf("DATA TYPE DOESN'T EXIST HOW DID YOU EVEN GET HERE?\n");
        }

        if(indicator_x ==1)
        {
          if(indicator_y==1)
          {
            Nconc_prospective[arm_i][i]++;
            Nloss_prospective[arm_i][i]++;

            Nconc_prospective[arm_i][n0+n1]++;
            Nwin_prospective[arm_i][n0+n1]++;
          }
          else if(indicator_y==0)
          {
            Ntied_conc_prospective[arm_i][i]++;
            Ntied_prospective[arm_i][i]++;

            Ntied_conc_prospective[arm_i][n0+n1]++;
            Ntied_prospective[arm_i][n0+n1]++;
          }
          else
          {
            Ndisc_prospective[arm_i][i]++;
            Nwin_prospective[arm_i][i]++;

            Ndisc_prospective[arm_i][n0+n1]++;
            Nloss_prospective[arm_i][n0+n1]++;
          }
        }
        else if(indicator_x == -1)
        {
          if(indicator_y==1)
          {
            Ndisc_prospective[arm_i][i]++;
            Nloss_prospective[arm_i][i]++;

            Ndisc_prospective[arm_i][n0+n1]++;
            Nwin_prospective[arm_i][n0+n1]++;
          }
          else if(indicator_y==0)
          {
            Ntied_conc_prospective[arm_i][i]++;
            Ntied_prospective[arm_i][i]++;

            Ntied_conc_prospective[arm_i][n0+n1]++;
            Ntied_prospective[arm_i][n0+n1]++;
          }
          else
          {
            Nconc_prospective[arm_i][i]++;
            Nwin_prospective[arm_i][i]++;

            Nconc_prospective[arm_i][n0+n1]++;
            Nloss_prospective[arm_i][n0+n1]++;
          }
        }
        else
        {
          if(indicator_y==1)
          {
            Nwin_prospective[arm_i][n0+n1]++;
            Nloss_prospective[arm_i][i]++;
          }
          else if(indicator_y==0)
          {
            Ntied_prospective[arm_i][n0+n1]++;
            Ntied_prospective[arm_i][i]++;
          }
          else
          {
            Nloss_prospective[arm_i][n0+n1]++;
            Nwin_prospective[arm_i][i]++;
          }
        }
      } // End for existing patients

      // The new patient ties with itself for the purposes of null
      Ntied_prospective[arm_i][n0+n1]++;
    } // End for prospective arm


    // Calculate genodds for both current and prospective imbalances

    std::vector<double> imbalance{0,0,0};

    genodds(imbalance,
            Nconc[cov_i][dummyCodeIndex],
            Ndisc[cov_i][dummyCodeIndex],
            Ntied_conc[cov_i][dummyCodeIndex],
            Nwin[cov_i][dummyCodeIndex],
            Nloss[cov_i][dummyCodeIndex],
            Ntied[cov_i][dummyCodeIndex],
            arm,
            n0,
            n1
    );
    out.imbalance.push_back(imbalance);

    std::vector<std::vector<double>> imbalance_prospective;

    for(int arm_i=0;arm_i<2;arm_i++)
    {
      std::vector<double> tmp_imbalance{0,0,0};

      genodds(tmp_imbalance,
              Nconc_prospective[arm_i],
              Ndisc_prospective[arm_i],
              Ntied_conc_prospective[arm_i],
              Nwin_prospective[arm_i],
              Nloss_prospective[arm_i],
              Ntied_prospective[arm_i],
              arm_prospective[arm_i],
              n0+(arm_i==0 ? 1 : 0),
              n1+(arm_i==0 ? 0 : 1)
      );

      // Don't make decision here, just store everything we'll need to make
      // a decision later

      imbalance_prospective.push_back(tmp_imbalance);

    }

    out.imbalance_prospective.push_back(imbalance_prospective);

  } // End for covariates

  //Imbalances are calculated, make intervention decision

  // Default is to not intervene
  out.biasDirection=-1;
  out.imbalanceStrength=0;

  if(voteMethod==0) // Simple majority
  {
    out.BD = -1;
    out.seBD = -1;

    // Count votes. If a vote is for control, de-increment
    // If a vote is for treatment, increment
    // majority direction it determined by sign of voteTally
    int voteTally = 0;

    // Considerations for weighted voting
    double voteStrength0=0;
    double voteStrength1=0;

    // Simple majority vote among significantly imbalanced covariates
    for(int cov_i=0; cov_i<nCovs; cov_i++)
    {
      // Z score for this covariate
      // For the sake of consistency, we're using SE rather than SEnull
      // this needs chased up later
      double thisTestStat = abs(log(out.imbalance[cov_i][0]) * out.imbalance[cov_i][0] /  out.imbalance[cov_i][1]);

      out.testStatistic[cov_i] = thisTestStat;

      if(thisTestStat > threshold[cov_i])
      {
        // Covariate is significant, check which direction we vote for
        if(abs(log(out.imbalance_prospective[cov_i][0][0])) < abs(log(out.imbalance_prospective[cov_i][1][0])))
        {
          // Arm 0 improves
          out.vote[cov_i] = 0;
          voteTally--;
          voteStrength0+= abs(log(out.imbalance[cov_i][0]));
        }
        else
        {
          // Arm 1 improves
          out.vote[cov_i] = 1;
          voteTally++;
          voteStrength1+= abs(log(out.imbalance[cov_i][0]));
        }


      }

    } // End for covariates

    if(voteTally>0)
    {
      out.biasDirection=1;
      out.imbalanceStrength=voteStrength1;
    }
    else if(voteTally<0)
    {
      out.biasDirection=0;
      out.imbalanceStrength=voteStrength0;
    }

  }
  else if(voteMethod==1) // Weighted majority
  {
    out.BD = -1;
    out.seBD = -1;

    // Weighted voting among significantly imbalanced covariates

    // As above, but now votes are weighted
    double voteTally = 0;

    // Considerations for weighted voting
    double voteStrength0=0;
    double voteStrength1=0;

    // Simple majority vote among significantly imbalanced covariates
    for(int cov_i=0; cov_i<nCovs; cov_i++)
    {
      // Z score for
      double thisTestStat = abs(log(out.imbalance[cov_i][0]) * out.imbalance[cov_i][0] /  out.imbalance[cov_i][1]);
      out.testStatistic[cov_i] = thisTestStat;

      if(thisTestStat > threshold[cov_i])
      {
        // Covariate is significant, check which direction we vote for
        if(abs(log(out.imbalance_prospective[cov_i][0][0])) < abs(log(out.imbalance_prospective[cov_i][1][0])))
        {
          // Arm 0 improves
          out.vote[cov_i] = 0;
          voteTally -= abs(log(out.imbalance[cov_i][0]));
          voteStrength0+= abs(log(out.imbalance[cov_i][0]));
        }
        else
        {
          // Arm 1 improves
          out.vote[cov_i] = 1;
          voteTally += abs(log(out.imbalance[cov_i][0]));
          voteStrength1+= abs(log(out.imbalance[cov_i][0]));
        }
      }

    } // End for covariate

    out.imbalanceStrength=abs(voteStrength1-voteStrength0);
    if(voteTally>0)
    {
      out.biasDirection=1;
    }
    else if(voteTally<0)
    {
      out.biasDirection=0;
    }


  }
  else if(voteMethod==2) // Pooled
  {
    // Pooled voting among all covariates

    double BD = 0;
    double seBD = 0;
    double BD_scale = 0;

    bool can_calculate = true;


    // Simple majority vote among significantly imbalanced covariates
    for(int cov_i=0; cov_i<nCovs; cov_i++)
    {

      // We only intervene when we can pool across all variables
      // This requires nonzero standard error
      // (also not getting {-1, -1, -1} error codes)
      if((out.imbalance[cov_i][1] <= 0) || (out.imbalance_prospective[cov_i][0][0] <= 0) || (out.imbalance_prospective[cov_i][1][0] <= 0) )
      {
        can_calculate = false;
        break;
      }

      // As written, we use SE, not SE_null for calculating BD.
      // This is probably something worth investigating in the future.

      double thisVariance = out.imbalance[cov_i][1] / out.imbalance[cov_i][0];

      thisVariance = thisVariance*thisVariance;

      // If assigning to control made things better, direction is negative
      int direction;
      if(abs(log(out.imbalance_prospective[cov_i][0][0])) < abs(log(out.imbalance_prospective[cov_i][1][0])))
      {
        direction=-1;
      }
      else if(abs(log(out.imbalance_prospective[cov_i][0][0])) > abs(log(out.imbalance_prospective[cov_i][1][0])))
      {
        direction=1;
      }
      else
      {
        direction=0;
      }

      out.vote[cov_i] = (direction==1 ? 1 : (direction==-1 ? 0 : -1));

      BD += weights[cov_i] * direction * abs(log(out.imbalance[cov_i][0])) / thisVariance;
      seBD += weights[cov_i] * weights[cov_i] / thisVariance;
      BD_scale += weights[cov_i]/thisVariance;

    }

    if(can_calculate)
    {

      BD = BD/BD_scale;
      seBD = sqrt(seBD/(BD_scale*BD_scale));

      out.BD = BD;
      out.seBD = seBD;

      double thisTestStat = abs(BD)/seBD;
      out.testStatistic[0] = thisTestStat;

      if(thisTestStat >= threshold[0])
      {
        out.imbalanceStrength=abs(BD);
        if(BD>0)
        {
          out.biasDirection=1;
        }
        else if(BD<0)
        {
          out.biasDirection=0;
        }
      }
    }
    else
    {

      for(int cov_i=0; cov_i<nCovs; cov_i++)
      {
        out.vote[cov_i]=-1;
      }

      out.BD = 0;
      out.seBD = -1;
      out.biasDirection=-1;
    }
  }
  else
  {
    Rprintf("Invalid voting method!\n");
  }

  // Intervention decision has been made!
  // Need to calculate randomisation probability.

  if(out.biasDirection==-1)
  {
    out.biasProbability_covariate=0.5;
  }
  else
  {
    out.biasProbability_covariate=getProbability(out.imbalanceStrength,biasMethod,parms);
  }

  // probability is standardised as probability of assignment to control group
  if(out.biasDirection==1) out.biasProbability_covariate = 1-out.biasProbability_covariate;


  out.biasProbability=out.biasProbability_covariate;



  // If group allocation imbalance is also being considered, evaluate this here.
  out.biasProbability_group=0.5;
  double groupImbalance = abs(((double) n0 - n1)/(n0+n1));


  if(groupImbalance > groupImbalanceThreshold){
    if(n0>n1)  out.biasProbability_group = 1-parms[0];
    if(n0<n1)  out.biasProbability_group = parms[0];
  }


  if(out.biasProbability_covariate!=out.biasProbability_group){

    if(out.biasProbability_covariate==0.5){ // covImbalance doesn't care

      out.biasProbability = out.biasProbability_group;

    } else {

      if(out.biasProbability_group==0.5){ // groupImbalance doesn't care

        out.biasProbability = out.biasProbability_covariate;

      } else {

        // # Both covariate balance and group imbalance care

        if(imbalancePriority  == 0){ // Abstain

          out.biasProbability = 0.5;

        } else if(imbalancePriority  == 1){  // Covariate

          out.biasProbability = out.biasProbability_covariate;

        } else if(imbalancePriority  == 2){  // Group

          out.biasProbability = out.biasProbability_group;

        }

      }
    }
  }


  return out;
}



/* * * * * * * * * * * * * * * * * * * *
 *
 *            EXPORTS TO R
 *
 * * * * * * * * * * * * * * * * * * * */

// [[Rcpp::export]]
List csmsb_get_decision(int nObs,
                  IntegerVector dataType_in,
                  IntegerVector nClasses_in,
                  NumericMatrix data,
                  IntegerVector arm_in,
                  NumericVector threshold_in,
                  NumericVector parms_in,
                  int voteMethod,
                  int biasMethod,
                  NumericVector weights_in,
                  double groupImbalanceThreshold,
                  int imbalancePriority
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

  if(voteMethod!=2){
    if(threshold_in.size()<nCov) Rprintf("Provide thresholds for all variables!\n");
  }

  std::vector<int> dataType(nCov);
  std::vector<int> nClasses(nCov);
  std::vector<double> threshold(nCov);
  std::vector<double> weights(nCov);
  std::vector<double> parms(parms_in.size());
  std::vector<std::vector<double>> covariates;
  std::vector<int> arm(nObs,0);
  int n0=0;
  int n1=0;

  std::vector<std::vector<std::vector<int>>> Nconc;
  std::vector<std::vector<std::vector<int>>> Ndisc;
  std::vector<std::vector<std::vector<int>>> Ntied_conc;
  std::vector<std::vector<std::vector<int>>> Nwin;
  std::vector<std::vector<std::vector<int>>> Nloss;
  std::vector<std::vector<std::vector<int>>> Ntied;

  // Initialise vectors
  for(int parm_i =0; parm_i < parms_in.size(); parm_i++){
    parms[parm_i] = parms_in[parm_i];
  }

  for(int cov_i=0;cov_i<nCov;cov_i++)
  {
    dataType[cov_i] = dataType_in[cov_i];
    nClasses[cov_i] = nClasses_in[cov_i];
    threshold[cov_i] = threshold_in[cov_i];
    weights[cov_i] = weights_in[cov_i];

    std::vector<double> thisColumn(nObs+1,-9999);
    covariates.push_back(thisColumn);

    std::vector<std::vector<int>> newNestedVector;
    for(int cov_j=0;cov_j<nClasses[cov_i];cov_j++)
    {
      // Nconc, etc will be initialised as a vector of vectors of vectors each of length nObs
      std::vector<int> emptyCount(nObs,0);
      newNestedVector.push_back(emptyCount);
    }

    Nconc.push_back(newNestedVector);
    Ndisc.push_back(newNestedVector);
    Ntied_conc.push_back(newNestedVector);
    Nwin.push_back(newNestedVector);
    Nloss.push_back(newNestedVector);
    Ntied.push_back(newNestedVector);
  }

  // Take existing data from R and feed it into our std::vector nonsense format.
  // Calculate concordance/discordance/etc at the same time
  // We go through this pain because it lets us keep track of concordant/discordant
  // pairs for genodds and cuts the complexity by a factor of O(n) when we run simulations


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

    // With the new patient added, update pair counts by comparing
    // with previous patients (also, test against itself!)
    for(int j=0;j<=i;j++)
    {
      for(int cov_i=0;cov_i<nCov;cov_i++)
      {
        // Need to iterate over all possible covariate values

        for(int class_i=0; class_i < nClasses[cov_i]; class_i++)
        {
          int indicator_x = arm[i] - arm[j];

          double epsilon = 0.0001;
          int indicator_y;
          if(dataType[cov_i]==0)
          {
            indicator_y = sgn(covariates[cov_i][i] - covariates[cov_i][j]);
          }
          else if(dataType[cov_i]==1)
          {

            indicator_y = (
              (abs(covariates[cov_i][i] - covariates[cov_i][j]) < epsilon) ?
            0 :
            sgn(covariates[cov_i][i] - covariates[cov_i][j])
            );
          }
          else if(dataType[cov_i]==2)
          {
            // Nominal variables get dummy coded
            // Check both i and j against class_i then treat as above
            int dummy_i = (abs(covariates[cov_i][i]-class_i) < epsilon) ? 1:0;
            int dummy_j = (abs(covariates[cov_i][j]-class_i) < epsilon) ? 1:0;
            indicator_y = dummy_i - dummy_j;
          }
          else
          {
            Rprintf("INVALID DATA TYPE\n");
          }

          if(indicator_x ==1)
          {
            if(indicator_y==1)
            {
              Nconc[cov_i][class_i][j]++;
              Nloss[cov_i][class_i][j]++;

              Nconc[cov_i][class_i][i]++;
              Nwin[cov_i][class_i][i]++;
            }
            else if(indicator_y==0)
            {
              Ntied_conc[cov_i][class_i][j]++;
              Ntied[cov_i][class_i][j]++;

              Ntied_conc[cov_i][class_i][i]++;
              Ntied[cov_i][class_i][i]++;
            }
            else
            {
              Ndisc[cov_i][class_i][j]++;
              Nwin[cov_i][class_i][j]++;

              Ndisc[cov_i][class_i][i]++;
              Nloss[cov_i][class_i][i]++;
            }
          }
          else if(indicator_x == -1)
          {
            if(indicator_y==1)
            {
              Ndisc[cov_i][class_i][j]++;
              Nloss[cov_i][class_i][j]++;

              Ndisc[cov_i][class_i][i]++;
              Nwin[cov_i][class_i][i]++;
            }
            else if(indicator_y==0)
            {
              Ntied_conc[cov_i][class_i][j]++;
              Ntied[cov_i][class_i][j]++;

              Ntied_conc[cov_i][class_i][i]++;
              Ntied[cov_i][class_i][i]++;
            }
            else
            {
              Nconc[cov_i][class_i][j]++;
              Nwin[cov_i][class_i][j]++;

              Nconc[cov_i][class_i][i]++;
              Nloss[cov_i][class_i][i]++;
            }
          }
          else
          {
            // We compare against ourself once
            if(i!=j)
            {
              if(indicator_y==1)
              {
                Nloss[cov_i][class_i][j]++;
                Nwin[cov_i][class_i][i]++;
              }
              else if(indicator_y==0)
              {
                Ntied[cov_i][class_i][j]++;
                Ntied[cov_i][class_i][i]++;
              }
              else
              {
                Nwin[cov_i][class_i][j]++;
                Nloss[cov_i][class_i][i]++;
              }
            }
            else
            {
              Ntied[cov_i][class_i][i]++;
            }
          }

        } // End class_i
      }  // end for cov_i
    } // end for j
  } // end for i


  // Do the same for the yet-to-be randomised patient
  // (don't update counts for concordance/discordance/etc)
  for(int cov_i=0;cov_i<nCov;cov_i++)
  {
    covariates[cov_i][n0+n1] = data(n0+n1,cov_i);
  }

  // Take all of this nonsense and get a decision
  csmsb_Intervention intervention = csmsb_interventionDecision(dataType,
                                                   nClasses,
                                                   covariates,
                                                   arm,
                                                   n0,
                                                   n1,
                                                   threshold,
                                                   parms,
                                                   voteMethod,
                                                   biasMethod,
                                                   weights,
                                                   groupImbalanceThreshold,
                                                   imbalancePriority,
                                                   Nconc,
                                                   Ndisc,
                                                   Ntied_conc,
                                                   Nwin,
                                                   Nloss,
                                                   Ntied
  );

  List out = List::create(
    Named("biasDirection") = intervention.biasDirection,
    Named("strength") = intervention.imbalanceStrength,
    Named("probability") = intervention.biasProbability,
    Named("probability_covariate") = intervention.biasProbability_covariate,
    Named("probability_group") = intervention.biasProbability_group,
    Named("imbalance") = wrap(intervention.imbalance),
    Named("imbalance_prospective") = wrap(intervention.imbalance_prospective),
    Named("vote") = wrap(intervention.vote),
    Named("testStatistic") = wrap(intervention.testStatistic),
    Named("BD")= intervention.BD,
    Named("seBD")=intervention.seBD
  );

  return out;
}


// [[Rcpp::export]]
List csmsb_simulate(IntegerVector dataType_in,
                    IntegerVector nClasses_in,
                    NumericMatrix data,
                    NumericVector roll_in,
                    NumericVector threshold_in,
                    NumericVector parms_in,
                    int voteMethod,
                    int biasMethod,
                    NumericVector weights_in,
                    double groupImbalanceThreshold,
                    int imbalancePriority
                    )
{

  int nObs = data.nrow();

  int nCov = dataType_in.length();

  std::vector<int> dataType(nCov);
  std::vector<int> nClasses(nCov);
  std::vector<double> threshold(nCov);
  std::vector<double> weights(nCov);
  std::vector<double> parms(parms_in.length());
  std::vector<double> roll(nObs);

  std::vector<std::vector<double>> covariates;
  std::vector<int> arm(nObs,0);

  std::vector<std::vector<std::vector<int>>> Nconc;
  std::vector<std::vector<std::vector<int>>> Ndisc;
  std::vector<std::vector<std::vector<int>>> Ntied_conc;
  std::vector<std::vector<std::vector<int>>> Nwin;
  std::vector<std::vector<std::vector<int>>> Nloss;
  std::vector<std::vector<std::vector<int>>> Ntied;

  for(int i=0;i<nObs;i++)
  {
    roll[i] = roll_in(i);
  }


  // Initialise vectors
  for(int parms_i=0;parms_i<parms_in.length();parms_i++)
  {
    parms[parms_i] = parms_in(parms_i);
  }

  for(int cov_i=0;cov_i<nCov;cov_i++)
  {
    dataType[cov_i] = dataType_in(cov_i);
    nClasses[cov_i] = nClasses_in(cov_i);
    threshold[cov_i] = threshold_in(cov_i);
    weights[cov_i] = weights_in[cov_i];

    std::vector<double> thisColumn(nObs+1,-9999);
    covariates.push_back(thisColumn);

    std::vector<std::vector<int>> newNestedVector;
    for(int cov_j=0;cov_j<nClasses[cov_i];cov_j++)
    {
      // Nconc, etc will be initialised as a vector of vectors of vectors each of length nObs
      std::vector<int> emptyCount(nObs,0);
      newNestedVector.push_back(emptyCount);
    }

    Nconc.push_back(newNestedVector);
    Ndisc.push_back(newNestedVector);
    Ntied_conc.push_back(newNestedVector);
    Nwin.push_back(newNestedVector);
    Nloss.push_back(newNestedVector);
    Ntied.push_back(newNestedVector);
  }

  // Randomisation outputs are defined here

  // Imbalance at each step for each covariate
  // SE of imbalance at each step for each covariate
  // Stored in paris of 2*i_cov and 2*i_cov+1
  NumericMatrix imbalance(nObs,2*nCov);

  // Vote placed by each covariate at each step
  IntegerMatrix vote(nObs,nCov);
  NumericMatrix testStatistic(nObs,nCov);

  // Randomisation probability (in favour of treatment)
  NumericVector randProb(nObs);
  NumericVector randProb_covariate(nObs);
  NumericVector randProb_group(nObs);

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

    // Randomise the new patient
    csmsb_Intervention intervention = csmsb_interventionDecision(dataType,
                                                     nClasses,
                                                     covariates,
                                                     arm,
                                                     n0,
                                                     n1,
                                                     threshold,
                                                     parms,
                                                     voteMethod,
                                                     biasMethod,
                                                     weights,
                                                     groupImbalanceThreshold,
                                                     imbalancePriority ,
                                                     Nconc,
                                                     Ndisc,
                                                     Ntied_conc,
                                                     Nwin,
                                                     Nloss,
                                                     Ntied
    );

    // Log randomisation and update current file

    for(int cov_i=0;cov_i<nCov;cov_i++)
    {
      imbalance(i,2*cov_i) = intervention.imbalance[cov_i][0];
      imbalance(i,2*cov_i+1) = intervention.imbalance[cov_i][1];
      vote(i,cov_i) = intervention.vote[cov_i];
      testStatistic(i,cov_i) = intervention.testStatistic[cov_i];
    }
    randProb(i) = intervention.biasProbability;
    randProb_covariate(i) = intervention.biasProbability_covariate;
    randProb_group(i) = intervention.biasProbability_group;

    int newArm;
    if(roll[i] < intervention.biasProbability)
    {
      newArm=0;
    }
    else
    {
      newArm=1;
    }

    arm[i] = newArm;
    if(arm[i]==0)
    {
      n0++;
    }
    else
    {
      n1++;
    }




    // With the new patient randomised, update pair counts by comparing
    // with previous patients (also, test against itself!)
    for(int j=0;j<=i;j++)
    {
      for(int cov_i=0;cov_i<nCov;cov_i++)
      {

        // Need to iterate over all possible covariate values

        for(int class_i=0; class_i < nClasses[cov_i]; class_i++)
        {


          int indicator_x = arm[i] - arm[j];

          double epsilon = 0.0001;
          int indicator_y;
          if(dataType[cov_i]==0)
          {
            indicator_y = sgn(covariates[cov_i][i] - covariates[cov_i][j]);
          }
          else if(dataType[cov_i]==1)
          {

            indicator_y = (
              (abs(covariates[cov_i][i] - covariates[cov_i][j]) < epsilon) ?
            0 :
            sgn(covariates[cov_i][i] - covariates[cov_i][j])
            );
          }
          else if(dataType[cov_i]==2)
          {
            // Nominal variables get dummy coded
            // Check both i and j against class_i then treat as above
            int dummy_i = (abs(covariates[cov_i][i]-class_i) < epsilon) ? 1:0;
            int dummy_j = (abs(covariates[cov_i][j]-class_i) < epsilon) ? 1:0;
            indicator_y = dummy_i - dummy_j;
          }
          else
          {
            Rprintf("INVALID DATA TYPE\n");
          }


          if(indicator_x ==1)
          {
            if(indicator_y==1)
            {
              Nconc[cov_i][class_i][j]++;
              Nloss[cov_i][class_i][j]++;

              Nconc[cov_i][class_i][i]++;
              Nwin[cov_i][class_i][i]++;
            }
            else if(indicator_y==0)
            {
              Ntied_conc[cov_i][class_i][j]++;
              Ntied[cov_i][class_i][j]++;

              Ntied_conc[cov_i][class_i][i]++;
              Ntied[cov_i][class_i][i]++;
            }
            else
            {
              Ndisc[cov_i][class_i][j]++;
              Nwin[cov_i][class_i][j]++;

              Ndisc[cov_i][class_i][i]++;
              Nloss[cov_i][class_i][i]++;
            }
          }
          else if(indicator_x == -1)
          {
            if(indicator_y==1)
            {
              Ndisc[cov_i][class_i][j]++;
              Nloss[cov_i][class_i][j]++;

              Ndisc[cov_i][class_i][i]++;
              Nwin[cov_i][class_i][i]++;
            }
            else if(indicator_y==0)
            {
              Ntied_conc[cov_i][class_i][j]++;
              Ntied[cov_i][class_i][j]++;

              Ntied_conc[cov_i][class_i][i]++;
              Ntied[cov_i][class_i][i]++;
            }
            else
            {
              Nconc[cov_i][class_i][j]++;
              Nwin[cov_i][class_i][j]++;

              Nconc[cov_i][class_i][i]++;
              Nloss[cov_i][class_i][i]++;
            }
          }
          else
          {
            // We compare against ourself once
            if(i!=j)
            {
              if(indicator_y==1)
              {
                Nloss[cov_i][class_i][j]++;
                Nwin[cov_i][class_i][i]++;
              }
              else if(indicator_y==0)
              {
                Ntied[cov_i][class_i][j]++;
                Ntied[cov_i][class_i][i]++;
              }
              else
              {
                Nwin[cov_i][class_i][j]++;
                Nloss[cov_i][class_i][i]++;
              }
            }
            else
            {
              Ntied[cov_i][class_i][i]++;
            }
          }

        } // End class_i
      }  // end for cov_i
    } // end for j



    // Add the next patient to the data frame should this not be the final randomisation
    if(i<(nObs-1))
    {
      for(int cov_i=0;cov_i<nCov;cov_i++)
      {
        covariates[cov_i][i+1] = data(i+1,cov_i);
      }
    }



  } // end for i



  end:

  List out = List::create(
    Named("arm")=arm,
    Named("imbalance")=imbalance,
    Named("vote")=vote,
    Named("testStatistic")=testStatistic,
    Named("randProb")=randProb,
    Named("randProb_covariate")=randProb_covariate,
    Named("randProb_group")=randProb_group
  );

  return out;
}

/* * * * * * * * * * * * * * * * * * * *
 *
 *          TEST FUNCTIONS
 *
 * * * * * * * * * * * * * * * * * * * */

// Test function exposing genodds

// [[Rcpp::export]]
NumericVector test_genodds(NumericVector Y,
                           IntegerVector X)
{

  // Run genodds on a given set of inputs
  // Should return the same values as the genodds package.

  // Get concordant, discordant, etc. pairings

  int nObs = X.length();

  std::vector<double> out(3,0);
  std::vector<int> Nconc(nObs,0);
  std::vector<int> Ndisc(nObs,0);
  std::vector<int> Ntied_conc(nObs,0);
  std::vector<int> Nwin(nObs,0);
  std::vector<int> Nloss(nObs,0);
  std::vector<int> Ntied(nObs,0);
  std::vector<int> arm(nObs,0);

  int n0=0;
  int n1=0;

  for(int i=0; i < nObs; i++)
  {

    arm[i] = X(i);
    if(arm[i]==0)
    {
      n0++;
    }
    else
    {
      n1++;
    }

    for(int j=0; j<nObs; j++)
    {
      // Copying code blocks from above is a cool and normal thing to do that
      // isn't going to cause bug-chasing issues at all

      int indicator_x = X(i) - X(j);
      int indicator_y = sgn(Y(i) - Y(j));

      if(indicator_x ==1)
      {
        if(indicator_y==1)
        {
          Nconc[i]++;
          Nwin[i]++;
        }
        else if(indicator_y==0)
        {
          Ntied_conc[i]++;
          Ntied[i]++;
        }
        else
        {
          Ndisc[i]++;
          Nloss[i]++;
        }
      }
      else if(indicator_x == -1)
      {
        if(indicator_y==1)
        {
          Ndisc[i]++;
          Nwin[i]++;
        }
        else if(indicator_y==0)
        {
          Ntied_conc[i]++;
          Ntied[i]++;
        }
        else
        {
          Nconc[i]++;
          Nloss[i]++;
        }
      }
      else
      {
        if(indicator_y==1)
        {
          Nwin[i]++;
        }
        else if(indicator_y==0)
        {
          Ntied[i]++;
        }
        else
        {
          Nloss[i]++;
        }
      }
    } // End for j
  } // End for i

//
//   for(int i=0; i<n0+n1;i++)
//   {
//     Rprintf("%d %d %d %d %d %d\n",
//             Nconc[i],
//             Ndisc[i],
//             Ntied_conc[i],
//             Nwin[i],
//             Nloss[i],
//             Ntied[i]
//     );
//   }

  genodds(out,
          Nconc,
          Ndisc,
          Ntied_conc,
          Nwin,
          Nloss,
          Ntied,
          arm,
          n0,
          n1
  );


  // In hindsight calling a temporary variable out was a bad idea
  NumericVector out_final_FINAL = Rcpp::wrap(out);

  return out_final_FINAL;

}

