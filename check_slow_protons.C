#include <TBranch.h>
#include <TCanvas.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TLine.h>
#include <TNamed.h>
#include <TObjArray.h>
#include <TPad.h>
#include <TParameter.h>
#include <TPaveText.h>
#include <TROOT.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <string>
#include <utility>
#include <vector>


enum class SupportClass {
  Unsupported = 0,
  Marginal = 1,
  Supported = 2
};


struct TimingShape {
  bool valid = false;
  bool boundHit = false;

  int fitStatus = -999;

  double kaonAmplitude = 0.0;
  double kaonAmplitudeError = 0.0;
  double kaonMean = 0.0;
  double kaonSigma = 0.0;

  double protonAmplitude = 0.0;
  double protonAmplitudeError = 0.0;
  double protonMean = 0.0;
  double protonSigma = 0.0;

  double otherAmplitude = 0.0;

  double separation = 0.0;
  double kaonSignificance = 0.0;
  double protonSignificance = 0.0;
  double poissonDeviance = 0.0;
  int goodnessNdf = 0;
  double poissonDevianceNdf = 0.0;
  double poissonDeviancePerEntry = 0.0;

  TF1 *fitFunction = nullptr;
};


struct SliceFitResult {
  bool valid = false;

  int fitStatus = -999;

  double kaonAmplitude = 0.0;
  double kaonAmplitudeError = 0.0;

  double protonAmplitude = 0.0;
  double protonAmplitudeError = 0.0;

  double otherAmplitude = 0.0;
  double otherAmplitudeError = 0.0;

  double kaonYield = 0.0;
  double protonYield = 0.0;
  double otherYield = 0.0;

  double dataYield = 0.0;
  double modelYield = 0.0;

  double modelDataRatio = 0.0;
  double poissonDeviance = 0.0;
  int goodnessNdf = 0;
  double poissonDevianceNdf = 0.0;
  double poissonDeviancePerEntry = 0.0;

  TF1 *fitFunction = nullptr;
};


struct PoissonGoodnessOfFit {
  double deviance = std::numeric_limits<double>::infinity();
  int ndf = 0;
  int fittedBins = 0;
  double fittedEntries = 0.0;
  double devianceNdf = std::numeric_limits<double>::infinity();
  double deviancePerEntry = std::numeric_limits<double>::infinity();
};


struct PeakPairSeed {
  bool valid = false;
  double lowerMean = 0.0;
  double upperMean = 0.0;
  double lowerHeight = 0.0;
  double upperHeight = 0.0;
};


std::string findAvailableBranch(
  TTree *tree,
  const std::vector<std::string> &candidates
) {
  if (!tree) {
    return "";
  }

  for (const std::string &candidate : candidates) {
    if (
      !candidate.empty() &&
      tree->GetBranch(candidate.c_str())
    ) {
      return candidate;
    }
  }

  return "";
}


void printRFLikeBranches(
  TTree *tree
) {
  if (!tree || !tree->GetListOfBranches()) {
    return;
  }

  bool printedHeader = false;

  TObjArray *branches = tree->GetListOfBranches();

  for (int index = 0; index < branches->GetEntries(); ++index) {
    auto *branch = dynamic_cast<TBranch *>(branches->At(index));

    if (!branch) {
      continue;
    }

    TString branchName = branch->GetName();
    TString lowerName = branchName;
    lowerName.ToLower();

    if (!lowerName.Contains("rf")) {
      continue;
    }

    if (!printedHeader) {
      std::cerr
        << "RF-like branches found in the tree:"
        << std::endl;
      printedHeader = true;
    }

    std::cerr
      << "  "
      << branchName
      << std::endl;
  }
}


const char *supportClassLabel(
  SupportClass supportClass
) {
  switch (supportClass) {
    case SupportClass::Supported:
      return "supported";

    case SupportClass::Marginal:
      return "marginal";

    default:
      return "unsupported";
  }
}


int supportClassColor(
  SupportClass supportClass
) {
  switch (supportClass) {
    case SupportClass::Supported:
      return kGreen + 2;

    case SupportClass::Marginal:
      return kOrange + 7;

    default:
      return kRed + 1;
  }
}


bool isNearBound(
  double value,
  double lowerBound,
  double upperBound,
  double fractionalTolerance
) {
  if (
    !std::isfinite(value) ||
    upperBound <= lowerBound
  ) {
    return true;
  }

  const double tolerance =
    fractionalTolerance * (upperBound - lowerBound);

  return
    value <= lowerBound + tolerance ||
    value >= upperBound - tolerance;
}


std::pair<double, double> findPeakSeed(
  const TH1D *histogram,
  double xMin,
  double xMax
) {
  if (!histogram) {
    return {
      0.0,
      0.5 * (xMin + xMax)
    };
  }

  const int firstBin = std::max(
    1,
    histogram->GetXaxis()->FindFixBin(xMin)
  );

  const int lastBin = std::min(
    histogram->GetNbinsX(),
    histogram->GetXaxis()->FindFixBin(
      std::nextafter(xMax, xMin)
    )
  );

  double maximum = 0.0;
  double maximumCenter =
    0.5 * (xMin + xMax);

  for (
    int bin = firstBin;
    bin <= lastBin;
    ++bin
  ) {
    const double content =
      histogram->GetBinContent(bin);

    if (content > maximum) {
      maximum = content;

      maximumCenter =
        histogram
          ->GetXaxis()
          ->GetBinCenter(bin);
    }
  }

  return {
    maximum,
    maximumCenter
  };
}


PeakPairSeed findSeparatedPeakPair(
  const TH1D *histogram,
  double xMin,
  double xMax,
  double minimumSeparation,
  double maximumSeparation
) {
  PeakPairSeed result;

  if (
    !histogram ||
    xMax <= xMin ||
    minimumSeparation <= 0.0 ||
    maximumSeparation <= minimumSeparation
  ) {
    return result;
  }

  const int firstBin = std::max(
    1,
    histogram->GetXaxis()->FindFixBin(xMin)
  );

  const int lastBin = std::min(
    histogram->GetNbinsX(),
    histogram->GetXaxis()->FindFixBin(
      std::nextafter(xMax, xMin)
    )
  );

  if (lastBin - firstBin < 4) {
    return result;
  }

  std::vector<double> smoothed(
    histogram->GetNbinsX() + 1,
    0.0
  );

  const int smoothingRadius = 2;

  for (int bin = firstBin; bin <= lastBin; ++bin) {
    double weightedSum = 0.0;
    double weightSum = 0.0;

    for (
      int offset = -smoothingRadius;
      offset <= smoothingRadius;
      ++offset
    ) {
      const int neighbor = bin + offset;

      if (neighbor < firstBin || neighbor > lastBin) {
        continue;
      }

      const double weight =
        static_cast<double>(smoothingRadius + 1 - std::abs(offset));

      weightedSum += weight * std::max(
        histogram->GetBinContent(neighbor),
        0.0
      );

      weightSum += weight;
    }

    if (weightSum > 0.0) {
      smoothed.at(bin) = weightedSum / weightSum;
    }
  }

  std::vector<int> peakBins;

  for (int bin = firstBin + 1; bin < lastBin; ++bin) {
    const double value = smoothed.at(bin);

    if (
      value > 0.0 &&
      value >= smoothed.at(bin - 1) &&
      value > smoothed.at(bin + 1)
    ) {
      peakBins.push_back(bin);
    }
  }

  // A very smooth or edge-truncated distribution can leave fewer than two
  // strict local maxima.  Add the strongest bins as fallback candidates;
  // the pair-separation requirement below still prevents one peak from
  // being counted twice.
  if (peakBins.size() < 2) {
    std::vector<int> allBins;

    for (int bin = firstBin; bin <= lastBin; ++bin) {
      allBins.push_back(bin);
    }

    std::sort(
      allBins.begin(),
      allBins.end(),
      [&smoothed](int left, int right) {
        return smoothed.at(left) > smoothed.at(right);
      }
    );

    const size_t maximumFallbackBins = std::min<size_t>(
      allBins.size(),
      24
    );

    for (size_t index = 0; index < maximumFallbackBins; ++index) {
      if (
        std::find(
          peakBins.begin(),
          peakBins.end(),
          allBins.at(index)
        ) == peakBins.end()
      ) {
        peakBins.push_back(allBins.at(index));
      }
    }
  }

  double bestScore = -1.0;

  for (size_t first = 0; first < peakBins.size(); ++first) {
    for (size_t second = first + 1; second < peakBins.size(); ++second) {
      const int firstPeakBin = peakBins.at(first);
      const int secondPeakBin = peakBins.at(second);

      const double firstMean =
        histogram->GetXaxis()->GetBinCenter(firstPeakBin);

      const double secondMean =
        histogram->GetXaxis()->GetBinCenter(secondPeakBin);

      const double separation = std::abs(secondMean - firstMean);

      if (
        separation < minimumSeparation ||
        separation > maximumSeparation
      ) {
        continue;
      }

      const double firstHeight = smoothed.at(firstPeakBin);
      const double secondHeight = smoothed.at(secondPeakBin);

      if (firstHeight <= 0.0 || secondHeight <= 0.0) {
        continue;
      }

      // Favor pairs for which both peaks are populated.  This suppresses
      // pairing a real peak with a small statistical ripple or with the
      // same particle peak in an adjacent RF period.
      const double score =
        std::min(firstHeight, secondHeight) *
        std::sqrt(firstHeight * secondHeight);

      if (score <= bestScore) {
        continue;
      }

      bestScore = score;
      result.valid = true;

      if (firstMean < secondMean) {
        result.lowerMean = firstMean;
        result.upperMean = secondMean;
        result.lowerHeight = firstHeight;
        result.upperHeight = secondHeight;
      } else {
        result.lowerMean = secondMean;
        result.upperMean = firstMean;
        result.lowerHeight = secondHeight;
        result.upperHeight = firstHeight;
      }
    }
  }

  return result;
}


std::pair<double, double> findCentralQuantileRange(
  const TH1D *histogram,
  double requestedMin,
  double requestedMax,
  double lowerTailFraction,
  double upperTailFraction,
  int paddingBins
) {
  if (
    !histogram ||
    requestedMax <= requestedMin ||
    lowerTailFraction < 0.0 ||
    upperTailFraction < 0.0 ||
    lowerTailFraction + upperTailFraction >= 1.0
  ) {
    return {requestedMin, requestedMax};
  }

  const int firstRequestedBin = std::max(
    1,
    histogram->GetXaxis()->FindFixBin(requestedMin)
  );

  const int lastRequestedBin = std::min(
    histogram->GetNbinsX(),
    histogram->GetXaxis()->FindFixBin(
      std::nextafter(requestedMax, requestedMin)
    )
  );

  const double total = histogram->Integral(
    firstRequestedBin,
    lastRequestedBin
  );

  if (total <= 0.0) {
    return {requestedMin, requestedMax};
  }

  const double lowerTarget = lowerTailFraction * total;
  const double upperTarget = (1.0 - upperTailFraction) * total;

  double cumulative = 0.0;
  int firstCentralBin = firstRequestedBin;
  int lastCentralBin = lastRequestedBin;
  bool foundLower = false;

  for (int bin = firstRequestedBin; bin <= lastRequestedBin; ++bin) {
    cumulative += std::max(histogram->GetBinContent(bin), 0.0);

    if (!foundLower && cumulative >= lowerTarget) {
      firstCentralBin = bin;
      foundLower = true;
    }

    if (cumulative >= upperTarget) {
      lastCentralBin = bin;
      break;
    }
  }

  firstCentralBin = std::max(
    firstRequestedBin,
    firstCentralBin - std::max(paddingBins, 0)
  );

  lastCentralBin = std::min(
    lastRequestedBin,
    lastCentralBin + std::max(paddingBins, 0)
  );

  const double rangeMin = std::max(
    requestedMin,
    histogram->GetXaxis()->GetBinLowEdge(firstCentralBin)
  );

  const double rangeMax = std::min(
    requestedMax,
    histogram->GetXaxis()->GetBinUpEdge(lastCentralBin)
  );

  if (rangeMax <= rangeMin) {
    return {requestedMin, requestedMax};
  }

  return {rangeMin, rangeMax};
}


std::pair<double, double> estimateTreeBranchCentralRange(
  TTree *tree,
  const std::string &branchName,
  const std::string &selectionCut,
  double fallbackMin,
  double fallbackMax
) {
  if (
    !tree ||
    branchName.empty() ||
    !tree->GetBranch(branchName.c_str()) ||
    fallbackMax <= fallbackMin
  ) {
    return {fallbackMin, fallbackMax};
  }

  tree->SetEstimate(
    std::max<Long64_t>(
      tree->GetEntries() + 1,
      1000000
    )
  );

  const Long64_t selectedRows = tree->Draw(
    branchName.c_str(),
    selectionCut.c_str(),
    "goff"
  );

  if (selectedRows < 25) {
    return {fallbackMin, fallbackMax};
  }

  const double *values = tree->GetV1();

  if (!values) {
    return {fallbackMin, fallbackMax};
  }

  std::vector<double> finiteValues;
  finiteValues.reserve(
    static_cast<size_t>(selectedRows)
  );

  for (Long64_t row = 0; row < selectedRows; ++row) {
    const double value = values[row];

    if (std::isfinite(value)) {
      finiteValues.push_back(value);
    }
  }

  if (finiteValues.size() < 25) {
    return {fallbackMin, fallbackMax};
  }

  std::sort(
    finiteValues.begin(),
    finiteValues.end()
  );

  auto quantileValue = [
    &finiteValues
  ](
    double quantile
  ) -> double {
    const double clampedQuantile = std::clamp(
      quantile,
      0.0,
      1.0
    );

    const size_t index = static_cast<size_t>(
      std::round(
        clampedQuantile *
        static_cast<double>(finiteValues.size() - 1)
      )
    );

    return finiteValues.at(
      std::min(
        index,
        finiteValues.size() - 1
      )
    );
  };

  double centralLow = quantileValue(0.001);
  double centralHigh = quantileValue(0.999);

  if (
    !std::isfinite(centralLow) ||
    !std::isfinite(centralHigh) ||
    centralHigh <= centralLow
  ) {
    return {fallbackMin, fallbackMax};
  }

  const double centralWidth = centralHigh - centralLow;
  const double padding = std::max(
    0.05 * centralWidth,
    1.0e-6
  );

  return {
    centralLow - padding,
    centralHigh + padding
  };
}


PoissonGoodnessOfFit computePoissonGoodnessOfFit(
  const TH1D *histogram,
  TF1 *fitFunction,
  double fitMin,
  double fitMax,
  int numberOfFreeParameters
) {
  PoissonGoodnessOfFit result;

  if (
    !histogram ||
    !fitFunction ||
    fitMax <= fitMin ||
    numberOfFreeParameters < 0
  ) {
    return result;
  }

  double deviance = 0.0;
  double fittedEntries = 0.0;
  int fittedBins = 0;

  for (
    int bin = 1;
    bin <= histogram->GetNbinsX();
    ++bin
  ) {
    const double binCenter =
      histogram->GetXaxis()->GetBinCenter(bin);

    if (
      binCenter < fitMin ||
      binCenter > fitMax
    ) {
      continue;
    }

    const double binLow =
      histogram->GetXaxis()->GetBinLowEdge(bin);

    const double binHigh =
      histogram->GetXaxis()->GetBinUpEdge(bin);

    const double binWidth = binHigh - binLow;

    if (!(binWidth > 0.0)) {
      continue;
    }

    const double observed = std::max(
      histogram->GetBinContent(bin),
      0.0
    );

    // The fit uses option I, so ROOT compares the bin content with the
    // bin-averaged function value rather than only the bin center.
    const double expected = std::max(
      fitFunction->Integral(binLow, binHigh) / binWidth,
      1.0e-12
    );

    if (observed > 0.0) {
      deviance += 2.0 * (
        expected -
        observed +
        observed * std::log(observed / expected)
      );
    } else {
      deviance += 2.0 * expected;
    }

    fittedEntries += observed;
    ++fittedBins;
  }

  result.deviance = deviance;
  result.fittedBins = fittedBins;
  result.fittedEntries = fittedEntries;
  result.ndf = fittedBins - numberOfFreeParameters;

  if (result.ndf > 0) {
    result.devianceNdf =
      result.deviance / result.ndf;
  }

  if (result.fittedEntries > 0.0) {
    result.deviancePerEntry =
      result.deviance / result.fittedEntries;
  }

  return result;
}


double evaluateGaussian(
  double x,
  double amplitude,
  double mean,
  double sigma
) {
  if (
    !std::isfinite(x) ||
    !std::isfinite(amplitude) ||
    !std::isfinite(mean) ||
    !std::isfinite(sigma) ||
    amplitude <= 0.0 ||
    sigma <= 0.0
  ) {
    return 0.0;
  }

  const double z =
    (x - mean) / sigma;

  return amplitude * std::exp(-0.5 * z * z);
}


double sumGaussianOverBins(
  const TH1D *histogram,
  double amplitude,
  double mean,
  double sigma,
  double fitMin,
  double fitMax
) {
  if (
    !histogram ||
    amplitude <= 0.0 ||
    sigma <= 0.0
  ) {
    return 0.0;
  }

  double yield = 0.0;

  for (
    int bin = 1;
    bin <= histogram->GetNbinsX();
    ++bin
  ) {
    const double x =
      histogram
        ->GetXaxis()
        ->GetBinCenter(bin);

    if (
      x < fitMin ||
      x > fitMax
    ) {
      continue;
    }

    yield += evaluateGaussian(
      x,
      amplitude,
      mean,
      sigma
    );
  }

  return yield;
}


double sumConstantOverBins(
  const TH1D *histogram,
  double amplitude,
  double fitMin,
  double fitMax
) {
  if (
    !histogram ||
    amplitude <= 0.0
  ) {
    return 0.0;
  }

  double yield = 0.0;

  for (
    int bin = 1;
    bin <= histogram->GetNbinsX();
    ++bin
  ) {
    const double x =
      histogram
        ->GetXaxis()
        ->GetBinCenter(bin);

    if (
      x < fitMin ||
      x > fitMax
    ) {
      continue;
    }

    yield += amplitude;
  }

  return yield;
}


int findEdgeBin(
  double value,
  const std::vector<double> &edges
) {
  if (
    edges.size() < 2 ||
    value < edges.front() ||
    value > edges.back()
  ) {
    return -1;
  }

  if (value == edges.back()) {
    return
      static_cast<int>(edges.size()) - 2;
  }

  const auto upper =
    std::upper_bound(
      edges.begin(),
      edges.end(),
      value
    );

  const int bin =
    static_cast<int>(
      std::distance(
        edges.begin(),
        upper
      )
    ) - 1;

  if (
    bin < 0 ||
    bin >= static_cast<int>(edges.size()) - 1
  ) {
    return -1;
  }

  return bin;
}


std::string makeDeltaCut(
  const std::string &deltaBranch,
  double low,
  double high,
  bool includeUpperEdge
) {
  if (includeUpperEdge) {
    return TString::Format(
      "%s >= %.17g && %s <= %.17g",
      deltaBranch.c_str(),
      low,
      deltaBranch.c_str(),
      high
    ).Data();
  }

  return TString::Format(
    "%s >= %.17g && %s < %.17g",
    deltaBranch.c_str(),
    low,
    deltaBranch.c_str(),
    high
  ).Data();
}


TimingShape fitGlobalTimingShape(
  TH1D *histogram,
  const std::string &functionName,
  double fitMin,
  double fitMax,
  double kaonMeanMin,
  double kaonMeanMax,
  double protonMeanMin,
  double protonMeanMax,
  bool protonPeakIsLower,
  double sigmaMin,
  double sigmaMax,
  double initialSigma,
  double minimumSeparation,
  double minimumAmplitudeSignificance,
  bool useDeviancePerEntryValidation,
  double maximumPoissonDevianceNdf,
  double maximumPoissonDeviancePerEntry,
  double boundFractionTolerance,
  int minimumEntries
) {
  TimingShape result;

  if (
    !histogram ||
    histogram->Integral() < minimumEntries
  ) {
    return result;
  }

  const auto kaonSeed =
    findPeakSeed(
      histogram,
      kaonMeanMin,
      kaonMeanMax
    );

  const auto protonSeed =
    findPeakSeed(
      histogram,
      protonMeanMin,
      protonMeanMax
    );

  const double histogramMaximum =
    std::max(
      histogram->GetMaximum(),
      1.0
    );

  auto *fitFunction = new TF1(
    functionName.c_str(),
    "[0] * exp(-0.5 * pow((x - [1]) / [2], 2))"
    " + "
    "[3] * exp(-0.5 * pow((x - [4]) / [5], 2))"
    " + [6]",
    fitMin,
    fitMax
  );

  fitFunction->SetParName(0, "K amplitude");
  fitFunction->SetParName(1, "K mean");
  fitFunction->SetParName(2, "K sigma");

  fitFunction->SetParName(3, "p amplitude");
  fitFunction->SetParName(4, "p mean");
  fitFunction->SetParName(5, "p sigma");

  fitFunction->SetParName(6, "other constant");

  fitFunction->SetParameter(
    0,
    std::max(
      kaonSeed.first,
      0.20 * histogramMaximum
    )
  );

  fitFunction->SetParLimits(
    0,
    0.0,
    100.0 * histogramMaximum
  );

  fitFunction->SetParameter(
    1,
    kaonSeed.second
  );

  fitFunction->SetParLimits(
    1,
    kaonMeanMin,
    kaonMeanMax
  );

  fitFunction->SetParameter(
    2,
    std::clamp(initialSigma, sigmaMin, sigmaMax)
  );

  fitFunction->SetParLimits(
    2,
    sigmaMin,
    sigmaMax
  );

  fitFunction->SetParameter(
    3,
    std::max(
      protonSeed.first,
      0.10 * histogramMaximum
    )
  );

  fitFunction->SetParLimits(
    3,
    0.0,
    100.0 * histogramMaximum
  );

  fitFunction->SetParameter(
    4,
    protonSeed.second
  );

  fitFunction->SetParLimits(
    4,
    protonMeanMin,
    protonMeanMax
  );

  fitFunction->SetParameter(
    5,
    std::clamp(initialSigma, sigmaMin, sigmaMax)
  );

  fitFunction->SetParLimits(
    5,
    sigmaMin,
    sigmaMax
  );

  fitFunction->SetParameter(
    6,
    0.02 * histogramMaximum
  );

  fitFunction->SetParLimits(
    6,
    0.0,
    10.0 * histogramMaximum
  );

  TFitResultPtr fitResult =
    histogram->Fit(
      fitFunction,
      "SRLIQ0"
    );

  result.fitStatus =
    static_cast<int>(fitResult);

  result.fitFunction =
    fitFunction;

  result.kaonAmplitude =
    fitFunction->GetParameter(0);

  result.kaonAmplitudeError =
    fitFunction->GetParError(0);

  result.kaonMean =
    fitFunction->GetParameter(1);

  result.kaonSigma =
    std::abs(
      fitFunction->GetParameter(2)
    );

  result.protonAmplitude =
    fitFunction->GetParameter(3);

  result.protonAmplitudeError =
    fitFunction->GetParError(3);

  result.protonMean =
    fitFunction->GetParameter(4);

  result.protonSigma =
    std::abs(
      fitFunction->GetParameter(5)
    );

  result.otherAmplitude =
    fitFunction->GetParameter(6);

  const PoissonGoodnessOfFit goodness =
    computePoissonGoodnessOfFit(
      histogram,
      fitFunction,
      fitMin,
      fitMax,
      7
    );

  result.poissonDeviance = goodness.deviance;
  result.goodnessNdf = goodness.ndf;
  result.poissonDevianceNdf = goodness.devianceNdf;
  result.poissonDeviancePerEntry =
    goodness.deviancePerEntry;

  const double separationDenominator =
    std::sqrt(
      result.kaonSigma * result.kaonSigma +
      result.protonSigma * result.protonSigma
    );

  if (separationDenominator > 0.0) {
    result.separation =
      std::abs(
        result.protonMean -
        result.kaonMean
      ) / separationDenominator;
  }

  if (result.kaonAmplitudeError > 0.0) {
    result.kaonSignificance =
      result.kaonAmplitude /
      result.kaonAmplitudeError;
  }

  if (result.protonAmplitudeError > 0.0) {
    result.protonSignificance =
      result.protonAmplitude /
      result.protonAmplitudeError;
  }

  result.boundHit =
    isNearBound(
      result.kaonMean,
      kaonMeanMin,
      kaonMeanMax,
      boundFractionTolerance
    ) ||
    isNearBound(
      result.protonMean,
      protonMeanMin,
      protonMeanMax,
      boundFractionTolerance
    ) ||
    isNearBound(
      result.kaonSigma,
      sigmaMin,
      sigmaMax,
      boundFractionTolerance
    ) ||
    isNearBound(
      result.protonSigma,
      sigmaMin,
      sigmaMax,
      boundFractionTolerance
    );

  result.valid =
    result.fitStatus == 0 &&
    !result.boundHit &&
    std::isfinite(result.kaonAmplitude) &&
    std::isfinite(result.protonAmplitude) &&
    std::isfinite(result.kaonMean) &&
    std::isfinite(result.protonMean) &&
    std::isfinite(result.kaonSigma) &&
    std::isfinite(result.protonSigma) &&
    std::isfinite(result.poissonDevianceNdf) &&
    std::isfinite(result.poissonDeviancePerEntry) &&
    result.kaonAmplitude > 0.0 &&
    result.protonAmplitude > 0.0 &&
    result.kaonSigma > 0.0 &&
    result.protonSigma > 0.0 &&
    (
      protonPeakIsLower
        ? result.protonMean < result.kaonMean
        : result.kaonMean < result.protonMean
    ) &&
    result.separation >= minimumSeparation &&
    result.kaonSignificance >=
      minimumAmplitudeSignificance &&
    result.protonSignificance >=
      minimumAmplitudeSignificance &&
    (
      useDeviancePerEntryValidation
        ? result.poissonDeviancePerEntry <=
            maximumPoissonDeviancePerEntry
        : result.poissonDevianceNdf <=
            maximumPoissonDevianceNdf
    );

  return result;
}


SliceFitResult fitDeltaTimingSlice(
  TH1D *histogram,
  const TimingShape &shape,
  const std::string &functionName,
  double fitMin,
  double fitMax,
  bool useDeviancePerEntryValidation,
  double maximumPoissonDevianceNdf,
  double maximumPoissonDeviancePerEntry,
  double minimumModelDataRatio,
  double maximumModelDataRatio,
  int minimumEntries
) {
  SliceFitResult result;

  if (
    !histogram ||
    !shape.valid ||
    histogram->Integral() < minimumEntries
  ) {
    return result;
  }

  const double histogramMaximum =
    std::max(
      histogram->GetMaximum(),
      1.0
    );

  const auto kaonSeed =
    findPeakSeed(
      histogram,
      shape.kaonMean -
        2.0 * shape.kaonSigma,
      shape.kaonMean +
        2.0 * shape.kaonSigma
    );

  const auto protonSeed =
    findPeakSeed(
      histogram,
      shape.protonMean -
        2.0 * shape.protonSigma,
      shape.protonMean +
        2.0 * shape.protonSigma
    );

  auto *fitFunction = new TF1(
    functionName.c_str(),
    "[0] * exp(-0.5 * pow((x - [1]) / [2], 2))"
    " + "
    "[3] * exp(-0.5 * pow((x - [4]) / [5], 2))"
    " + [6]",
    fitMin,
    fitMax
  );

  fitFunction->SetParName(0, "K amplitude");
  fitFunction->SetParName(1, "K mean");
  fitFunction->SetParName(2, "K sigma");

  fitFunction->SetParName(3, "p amplitude");
  fitFunction->SetParName(4, "p mean");
  fitFunction->SetParName(5, "p sigma");

  fitFunction->SetParName(6, "other constant");

  fitFunction->SetParameter(
    0,
    std::max(
      kaonSeed.first,
      0.10 * histogramMaximum
    )
  );

  fitFunction->SetParLimits(
    0,
    0.0,
    100.0 * histogramMaximum
  );

  fitFunction->FixParameter(
    1,
    shape.kaonMean
  );

  fitFunction->FixParameter(
    2,
    shape.kaonSigma
  );

  fitFunction->SetParameter(
    3,
    std::max(
      protonSeed.first,
      0.05 * histogramMaximum
    )
  );

  fitFunction->SetParLimits(
    3,
    0.0,
    100.0 * histogramMaximum
  );

  fitFunction->FixParameter(
    4,
    shape.protonMean
  );

  fitFunction->FixParameter(
    5,
    shape.protonSigma
  );

  fitFunction->SetParameter(
    6,
    0.02 * histogramMaximum
  );

  fitFunction->SetParLimits(
    6,
    0.0,
    10.0 * histogramMaximum
  );

  TFitResultPtr fitResult =
    histogram->Fit(
      fitFunction,
      "SRLIQ0"
    );

  result.fitStatus =
    static_cast<int>(fitResult);

  result.fitFunction =
    fitFunction;

  result.kaonAmplitude =
    fitFunction->GetParameter(0);

  result.kaonAmplitudeError =
    fitFunction->GetParError(0);

  result.protonAmplitude =
    fitFunction->GetParameter(3);

  result.protonAmplitudeError =
    fitFunction->GetParError(3);

  result.otherAmplitude =
    fitFunction->GetParameter(6);

  result.otherAmplitudeError =
    fitFunction->GetParError(6);

  result.kaonYield =
    sumGaussianOverBins(
      histogram,
      result.kaonAmplitude,
      shape.kaonMean,
      shape.kaonSigma,
      fitMin,
      fitMax
    );

  result.protonYield =
    sumGaussianOverBins(
      histogram,
      result.protonAmplitude,
      shape.protonMean,
      shape.protonSigma,
      fitMin,
      fitMax
    );

  result.otherYield =
    sumConstantOverBins(
      histogram,
      result.otherAmplitude,
      fitMin,
      fitMax
    );

  const int firstFitBin =
    std::max(
      1,
      histogram
        ->GetXaxis()
        ->FindFixBin(fitMin)
    );

  const int lastFitBin =
    std::min(
      histogram->GetNbinsX(),
      histogram
        ->GetXaxis()
        ->FindFixBin(
          std::nextafter(
            fitMax,
            fitMin
          )
        )
    );

  result.dataYield =
    histogram->Integral(
      firstFitBin,
      lastFitBin
    );

  result.modelYield =
    result.kaonYield +
    result.protonYield +
    result.otherYield;

  if (result.dataYield > 0.0) {
    result.modelDataRatio =
      result.modelYield /
      result.dataYield;
  }

  const PoissonGoodnessOfFit goodness =
    computePoissonGoodnessOfFit(
      histogram,
      fitFunction,
      fitMin,
      fitMax,
      3
    );

  result.poissonDeviance = goodness.deviance;
  result.goodnessNdf = goodness.ndf;
  result.poissonDevianceNdf = goodness.devianceNdf;
  result.poissonDeviancePerEntry =
    goodness.deviancePerEntry;

  result.valid =
    result.fitStatus == 0 &&
    std::isfinite(result.kaonAmplitude) &&
    std::isfinite(result.protonAmplitude) &&
    std::isfinite(result.otherAmplitude) &&
    std::isfinite(result.kaonYield) &&
    std::isfinite(result.protonYield) &&
    std::isfinite(result.otherYield) &&
    std::isfinite(result.modelDataRatio) &&
    std::isfinite(result.poissonDevianceNdf) &&
    std::isfinite(result.poissonDeviancePerEntry) &&
    result.kaonAmplitude >= 0.0 &&
    result.protonAmplitude >= 0.0 &&
    result.otherAmplitude >= 0.0 &&
    result.modelYield > 0.0 &&
    result.modelDataRatio >=
      minimumModelDataRatio &&
    result.modelDataRatio <=
      maximumModelDataRatio &&
    (
      useDeviancePerEntryValidation
        ? result.poissonDeviancePerEntry <=
            maximumPoissonDeviancePerEntry
        : result.poissonDevianceNdf <=
            maximumPoissonDevianceNdf
    );

  return result;
}


double evaluateEventProtonProbability(
  double timing,
  const TimingShape &shape,
  const SliceFitResult &fitResult,
  double denominatorFloor
) {
  if (
    !shape.valid ||
    !fitResult.valid
  ) {
    return 0.0;
  }

  const double protonValue =
    evaluateGaussian(
      timing,
      fitResult.protonAmplitude,
      shape.protonMean,
      shape.protonSigma
    );

  const double kaonValue =
    evaluateGaussian(
      timing,
      fitResult.kaonAmplitude,
      shape.kaonMean,
      shape.kaonSigma
    );

  const double otherValue =
    std::max(
      fitResult.otherAmplitude,
      0.0
    );

  const double denominator =
    protonValue +
    kaonValue +
    otherValue;

  if (
    !std::isfinite(denominator) ||
    denominator <= denominatorFloor
  ) {
    return 0.0;
  }

  return std::clamp(
    protonValue / denominator,
    0.0,
    1.0
  );
}


TF1 *makeGaussianComponent(
  const std::string &name,
  double amplitude,
  double mean,
  double sigma,
  double xMin,
  double xMax,
  int color
) {
  auto *function = new TF1(
    name.c_str(),
    "[0] * exp(-0.5 * pow((x - [1]) / [2], 2))",
    xMin,
    xMax
  );

  function->SetParameter(0, amplitude);
  function->SetParameter(1, mean);
  function->SetParameter(2, sigma);

  function->SetLineColor(color);
  function->SetLineWidth(2);
  function->SetLineStyle(2);

  return function;
}


TF1 *makeConstantComponent(
  const std::string &name,
  double amplitude,
  double xMin,
  double xMax,
  int color
) {
  auto *function = new TF1(
    name.c_str(),
    "[0]",
    xMin,
    xMax
  );

  function->SetParameter(
    0,
    std::max(
      amplitude,
      0.0
    )
  );

  function->SetLineColor(color);
  function->SetLineWidth(2);
  function->SetLineStyle(3);

  return function;
}


void check_slow_protons(
  const char *phi_setting = "Left",
  const char *Q2 = "3p0",
  const char *W = "3p14",
  const char *eps_setting = "low"
) {
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  const char *macroVersion = "check_slow_protons.9";

  std::cout
    << "Running "
    << macroVersion
    << std::endl;

  const TString phiSetting =
    phi_setting ? phi_setting : "";

  const TString q2Setting =
    Q2 ? Q2 : "";

  const TString wSetting =
    W ? W : "";

  const TString epsSetting =
    eps_setting ? eps_setting : "";

  if (
    phiSetting.IsNull() ||
    q2Setting.IsNull() ||
    wSetting.IsNull() ||
    epsSetting.IsNull()
  ) {
    std::cerr
      << "Usage: check_slow_protons("
      << "\"<phi>\", \"<Q2>\", "
      << "\"<W>\", \"<eps>\")"
      << std::endl;

    return;
  }

  const TString kinematicSetting =
    TString::Format(
      "Q%sW%s",
      q2Setting.Data(),
      wSetting.Data()
    );

  const TString inputFileName =
    TString::Format(
      "%s_kaon_Analysed_Data_Q%sW%s_%se.root",
      phiSetting.Data(),
      q2Setting.Data(),
      wSetting.Data(),
      epsSetting.Data()
    );

  const TString inputDirectory =
    TString::Format(
      "$cache_transfer/%s/Trial_49/"
      "2026June29_H00M23S53/root",
      kinematicSetting.Data()
    );

  TString filename = gSystem->ExpandPathName(
    TString::Format(
      "%s/%s",
      inputDirectory.Data(),
      inputFileName.Data()
    )
  );

  const std::string treeName =
    "Cut_Kaon_Events_prompt_noRF";

  const std::string outputBase =
    TString::Format(
      "%s_kaon_proton_cleaning_Q%sW%s_%se",
      phiSetting.Data(),
      q2Setting.Data(),
      wSetting.Data(),
      epsSetting.Data()
    ).Data();

  const std::string outputPDF =
    outputBase + ".pdf";

  const std::string outputROOT =
    outputBase + ".root";

  TFile *inputFile =
    TFile::Open(
      filename,
      "READ"
    );

  if (
    !inputFile ||
    inputFile->IsZombie()
  ) {
    std::cerr
      << "Unable to open input file: "
      << filename
      << std::endl;

    return;
  }

  TTree *tree =
    inputFile->Get<TTree>(
      treeName.c_str()
    );

  if (!tree) {
    std::cerr
      << "Tree "
      << treeName
      << " not found."
      << std::endl;

    inputFile->ls();
    inputFile->Close();
    return;
  }

  const std::string aeroBranch =
    "P_aero_npeSum";

  const std::string ctTimeBranch =
    "CTime_ROC1";

  std::string timeBranch =
    ctTimeBranch;

  std::string rfTimeBranch;
  bool rfTimingAttempted = false;
  bool rfTimingSelected = false;
  bool ctTimingEvaluated = false;
  bool ctFallbackUsed = false;  // retained as compatibility metadata; v9 does not use fallback selection
  int rfProbeValidShapes = 0;
  int ctProbeValidShapes = 0;
  std::string timingSelectionReason;

  const std::string deltaBranch =
    "ssdelta";

  const std::string mmBranch =
    "MM";

  const std::vector<std::string> requiredBranches = {
    aeroBranch,
    ctTimeBranch,
    deltaBranch,
    mmBranch,
    "ssxptar",
    "ssyptar",
    "hsdelta",
    "hsxptar",
    "hsyptar"
  };

  for (
    const std::string &branch :
    requiredBranches
  ) {
    if (!tree->GetBranch(branch.c_str())) {
      std::cerr
        << "Required branch not found: "
        << branch
        << std::endl;

      inputFile->Close();
      return;
    }
  }

  // ------------------------------------------------------------------
  // Configuration.
  // ------------------------------------------------------------------

  const double aeroMin = 0.0;
  const double aeroMax = 25.0;

  const bool isHighEpsilon =
    epsSetting.EqualTo(
      "high",
      TString::kIgnoreCase
    );

  // The input tree remains the uncut *_noRF tree.  Every normal run
  // evaluates both the best available RF timing branch and CTime_ROC1,
  // then chooses the timing variable with the stronger validated global
  // proton-kaon discrimination.  No RF phase correction or wrapping is
  // applied.
  const double beamBunchSpacingNs =
    isHighEpsilon ? 4.0 : 2.0;

  // CT uses the same procedure for low and high epsilon.  High epsilon
  // only retains the wider display range implied by the 4 ns bunch spacing.
  const double ctTimeMin = -beamBunchSpacingNs;
  const double ctTimeMax = beamBunchSpacingNs;
  const double ctTimingFitMin = ctTimeMin;
  const double ctTimingFitMax = ctTimeMax;

  const double ctKaonMeanMin = -0.45;
  const double ctKaonMeanMax = 0.20;
  const double ctProtonMeanMin = 0.20;
  const double ctProtonMeanMax = 0.95;

  double timeMin = ctTimeMin;
  double timeMax = ctTimeMax;

  double timingFitMin = ctTimingFitMin;
  double timingFitMax = ctTimingFitMax;

  double kaonMeanMin = ctKaonMeanMin;
  double kaonMeanMax = ctKaonMeanMax;
  double protonMeanMin = ctProtonMeanMin;
  double protonMeanMax = ctProtonMeanMax;

  const double timingSigmaMin = 0.03;
  double timingSigmaMax = 0.45;
  double timingSigmaInitial = 0.15;

  const double minimumGlobalSeparation = 0.75;
  const double minimumGlobalAmplitudeSignificance = 2.0;

  // Use the same reduced Poisson-deviance validation for both epsilon
  // settings.
  const bool useDeviancePerEntryValidation = false;

  const double maximumGlobalPoissonDevianceNdf = 5.0;
  const double maximumGlobalPoissonDeviancePerEntry = 0.85;
  const double globalBoundFractionTolerance = 0.02;

  const double maximumSlicePoissonDevianceNdf = 5.0;
  const double maximumSlicePoissonDeviancePerEntry = 1.00;
  const double minimumSliceModelDataRatio = 0.50;
  const double maximumSliceModelDataRatio = 1.50;

  const double deltaMin = -10.0;
  const double deltaMax = 20.0;

  const int nDeltaBins = 10;

  const double mmMin = 0.70;
  const double mmMax = 1.50;

  const int nMMBins = 160;

  const int nAeroHistogramBins = 75;
  // Preserve approximately the original 0.0305 ns timing-bin width.
  // Eight ns gives 262 bins; four ns gives 131 bins.
  const int ctTimeHistogramBins =
    isHighEpsilon ? 262 : 131;

  int nTimeHistogramBins =
    ctTimeHistogramBins;

  const std::vector<double> aeroEdges = {
    0.0,
    3.0,
    6.0,
    10.0,
    15.0,
    25.0
  };

  const int nAeroSlices =
    static_cast<int>(
      aeroEdges.size()
    ) - 1;

  const int minimumGlobalSliceEntries = 200;
  const int minimumDeltaSliceEntries = 30;

  const int minimumSupportedSlices = 2;
  const int minimumMarginalSlices = 1;

  const double minimumSupportedCoverage = 0.35;
  const double minimumMarginalCoverage = 0.15;

  const double minimumModeledYield = 5.0;

  const double eventProbabilityDenominatorFloor =
    1.0e-12;

  const double marginalUncertaintyInflation = 2.0;

  const double lowMMMin = 0.80;
  const double lowMMMax = 0.90;

  const double lambdaMMMin = 1.105;
  const double lambdaMMMax = 1.125;

  const double deltaBinWidth =
    (deltaMax - deltaMin) /
    nDeltaBins;

  // ------------------------------------------------------------------
  // Analysis cuts and RF/CT timing comparison.
  // ------------------------------------------------------------------

  const std::string acceptanceCut =
    "ssdelta >= -10.0 && ssdelta <= 20.0 && "
    "ssxptar >= -0.06 && ssxptar <= 0.06 && "
    "ssyptar >= -0.04 && ssyptar <= 0.04 && "
    "hsdelta >= -8.0 && hsdelta <= 8.0 && "
    "hsxptar >= -0.08 && hsxptar <= 0.08 && "
    "hsyptar >= -0.045 && hsyptar <= 0.045";

  const char *disableRFEnvironment =
    gSystem->Getenv("PROTON_CHECKER_DISABLE_RF");

  const bool disableRFTiming =
    disableRFEnvironment &&
    std::string(disableRFEnvironment) == "1";

  std::vector<std::string> rfBranchCandidates;

  const char *forcedRFBranch =
    gSystem->Getenv("PROTON_CHECKER_RF_BRANCH");

  if (
    !disableRFTiming &&
    forcedRFBranch &&
    forcedRFBranch[0] != '\0'
  ) {
    rfBranchCandidates.emplace_back(forcedRFBranch);
  }

  const std::vector<std::string> defaultRFBranchCandidates = {
    "RF",
    "RFTime",
    "RF_time",
    "RFTime_ROC1",
    "P_RFTime",
    "P_RF_tdcTime",
    "P_RF_adcTime",
    "RF_Dist",
    "RF_Distance",
    "P_RF_Dist",
    "P_RF_Distance",
    "P_RF_Dist_Track",
    "H_RFTime",
    "H_RF_tdcTime"
  };

  for (
    const std::string &candidate :
    defaultRFBranchCandidates
  ) {
    if (disableRFTiming) {
      break;
    }
    if (
      std::find(
        rfBranchCandidates.begin(),
        rfBranchCandidates.end(),
        candidate
      ) == rfBranchCandidates.end()
    ) {
      rfBranchCandidates.push_back(candidate);
    }
  }

  struct TimingBranchProbe {
    std::string branch;
    int validShapes = 0;
    double displayMin = 0.0;
    double displayMax = 0.0;
    double fitMin = 0.0;
    double fitMax = 0.0;
    double kaonMeanMin = 0.0;
    double kaonMeanMax = 0.0;
    double protonMeanMin = 0.0;
    double protonMeanMax = 0.0;
    double sigmaMax = 0.0;
    double sigmaInitial = 0.0;
    int histogramBins = 0;
    bool peakPairFound = false;
    double protonSeedMean = 0.0;
    double kaonSeedMean = 0.0;
    double meanSeparation = 0.0;
    double meanPoissonDevianceNdf =
      std::numeric_limits<double>::infinity();
  };

  auto countValidGlobalShapesForBranch = [
    &
  ](
    const std::string &candidateBranch,
    TimingBranchProbe &probe
  ) -> int {
    if (
      candidateBranch.empty() ||
      !tree->GetBranch(candidateBranch.c_str())
    ) {
      return 0;
    }

    const std::string candidateBaseCut =
      TString::Format(
        "%s >= %.17g && %s <= %.17g",
        aeroBranch.c_str(),
        aeroMin,
        aeroBranch.c_str(),
        aeroMax
      ).Data();

    const std::string candidatePreselectionCut =
      "(" + acceptanceCut + ") && (" +
      candidateBaseCut + ")";

    const auto candidateDisplayRange =
      estimateTreeBranchCentralRange(
        tree,
        candidateBranch,
        candidatePreselectionCut,
        -beamBunchSpacingNs,
        beamBunchSpacingNs
      );

    probe.branch = candidateBranch;
    probe.displayMin = candidateDisplayRange.first;
    probe.displayMax = candidateDisplayRange.second;
    probe.fitMin = probe.displayMin;
    probe.fitMax = probe.displayMax;

    if (probe.displayMax <= probe.displayMin) {
      return 0;
    }

    // Search the unmodified RF variable for a proton-kaon peak pair.  For
    // low epsilon, P_RF_Dist can span roughly two 2-ns RF periods.  Fitting
    // the entire 4-ns span duplicates the particle peaks and makes the
    // two-Gaussian model ill-defined.  Select one contiguous bunch-period
    // window around the strongest two-peak structure; no RF value is
    // shifted, corrected, or wrapped.
    const double rawDisplayMin = probe.displayMin;
    const double rawDisplayMax = probe.displayMax;
    const double rawDisplayWidth = rawDisplayMax - rawDisplayMin;

    const int searchBins = std::max(
      160,
      static_cast<int>(
        std::round(rawDisplayWidth / 0.015)
      )
    );

    gROOT->cd();

    auto *rfWindowSearch = new TH1D(
      TString::Format(
        "h_rf_window_search_%s",
        candidateBranch.c_str()
      ),
      "RF window search",
      searchBins,
      rawDisplayMin,
      rawDisplayMax
    );

    rfWindowSearch->Sumw2();

    tree->Draw(
      TString::Format(
        "%s>>%s",
        candidateBranch.c_str(),
        rfWindowSearch->GetName()
      ),
      candidatePreselectionCut.c_str(),
      "goff"
    );

    rfWindowSearch->SetDirectory(nullptr);

    const PeakPairSeed rawPeakPair = findSeparatedPeakPair(
      rfWindowSearch,
      rawDisplayMin,
      rawDisplayMax,
      std::max(0.18, 0.10 * beamBunchSpacingNs),
      0.80 * beamBunchSpacingNs
    );

    if (
      rawPeakPair.valid &&
      rawDisplayWidth > 1.15 * beamBunchSpacingNs
    ) {
      const double pairCenter =
        0.5 * (rawPeakPair.lowerMean + rawPeakPair.upperMean);

      const double selectedWidth = beamBunchSpacingNs;
      double selectedMin = pairCenter - 0.5 * selectedWidth;
      double selectedMax = pairCenter + 0.5 * selectedWidth;

      if (selectedMin < rawDisplayMin) {
        selectedMax += rawDisplayMin - selectedMin;
        selectedMin = rawDisplayMin;
      }

      if (selectedMax > rawDisplayMax) {
        selectedMin -= selectedMax - rawDisplayMax;
        selectedMax = rawDisplayMax;
      }

      probe.displayMin = std::max(selectedMin, rawDisplayMin);
      probe.displayMax = std::min(selectedMax, rawDisplayMax);
    }

    delete rfWindowSearch;

    const double candidateWidth =
      probe.displayMax - probe.displayMin;

    const int candidateBins = std::max(
      80,
      static_cast<int>(
        std::round(
          candidateWidth / 0.0305
        )
      )
    );

    probe.histogramBins = candidateBins;

    const std::string candidateRangeCut =
      TString::Format(
        "%s >= %.17g && %s <= %.17g && "
        "%s >= %.17g && %s <= %.17g",
        aeroBranch.c_str(),
        aeroMin,
        aeroBranch.c_str(),
        aeroMax,
        candidateBranch.c_str(),
        probe.displayMin,
        candidateBranch.c_str(),
        probe.displayMax
      ).Data();

    const std::string candidateAnalysisCut =
      "(" + acceptanceCut + ") && (" +
      candidateRangeCut + ")";

    gROOT->cd();

    auto *probePID = new TH2D(
      TString::Format(
        "h_rf_first_probe_pid_%s",
        candidateBranch.c_str()
      ),
      "RF global timing probe",
      nAeroHistogramBins,
      aeroMin,
      aeroMax,
      candidateBins,
      probe.displayMin,
      probe.displayMax
    );

    probePID->Sumw2();

    tree->Draw(
      TString::Format(
        "%s:%s>>%s",
        candidateBranch.c_str(),
        aeroBranch.c_str(),
        probePID->GetName()
      ),
      candidateAnalysisCut.c_str(),
      "goff"
    );

    probePID->SetDirectory(nullptr);

    if (probePID->Integral() <= 0.0) {
      delete probePID;
      return 0;
    }

    auto *allAeroTiming = probePID->ProjectionY(
      TString::Format(
        "h_rf_first_probe_all_aero_%s",
        candidateBranch.c_str()
      ),
      1,
      probePID->GetNbinsX(),
      "e"
    );

    allAeroTiming->SetDirectory(nullptr);

    const auto centralFitRange = findCentralQuantileRange(
      allAeroTiming,
      probe.displayMin,
      probe.displayMax,
      5.0e-4,
      5.0e-4,
      2
    );

    probe.fitMin = centralFitRange.first;
    probe.fitMax = centralFitRange.second;

    if (probe.fitMax <= probe.fitMin) {
      delete allAeroTiming;
      delete probePID;
      return 0;
    }

    const PeakPairSeed fittedPeakPair = findSeparatedPeakPair(
      allAeroTiming,
      probe.fitMin,
      probe.fitMax,
      std::max(0.18, 0.10 * beamBunchSpacingNs),
      0.80 * beamBunchSpacingNs
    );

    double split = 0.5 * (probe.fitMin + probe.fitMax);

    if (fittedPeakPair.valid) {
      split = 0.5 * (
        fittedPeakPair.lowerMean +
        fittedPeakPair.upperMean
      );

      probe.peakPairFound = true;
      probe.protonSeedMean = fittedPeakPair.lowerMean;
      probe.kaonSeedMean = fittedPeakPair.upperMean;
    }

    delete allAeroTiming;

    // RF timing has the proton peak on the lower-time side and the kaon
    // peak on the upper-time side.  Split the allowed means between the
    // two data-driven peak seeds rather than at the center of the raw RF
    // range.  Use the same sigma constraints and seeds as the CT fit.
    probe.protonMeanMin = probe.fitMin;
    probe.protonMeanMax = split;
    probe.kaonMeanMin = split;
    probe.kaonMeanMax = probe.fitMax;
    probe.sigmaMax = timingSigmaMax;
    probe.sigmaInitial = timingSigmaInitial;

    int validShapes = 0;
    double validSeparationSum = 0.0;
    double validDevianceNdfSum = 0.0;

    for (
      int aeroSlice = 0;
      aeroSlice < nAeroSlices;
      ++aeroSlice
    ) {
      const double aeroLow =
        aeroEdges.at(aeroSlice);

      const double aeroHigh =
        aeroEdges.at(aeroSlice + 1);

      const int firstXBin = std::max(
        1,
        probePID
          ->GetXaxis()
          ->FindFixBin(
            std::nextafter(
              aeroLow,
              aeroHigh
            )
          )
      );

      const int lastXBin = std::min(
        probePID->GetNbinsX(),
        probePID
          ->GetXaxis()
          ->FindFixBin(
            std::nextafter(
              aeroHigh,
              aeroLow
            )
          )
      );

      auto *projection = probePID->ProjectionY(
        TString::Format(
          "h_rf_first_probe_aero_%s_%d",
          candidateBranch.c_str(),
          aeroSlice
        ),
        firstXBin,
        lastXBin,
        "e"
      );

      projection->SetDirectory(nullptr);

      TimingShape shape = fitGlobalTimingShape(
        projection,
        TString::Format(
          "f_rf_first_probe_aero_%s_%d",
          candidateBranch.c_str(),
          aeroSlice
        ).Data(),
        probe.fitMin,
        probe.fitMax,
        probe.kaonMeanMin,
        probe.kaonMeanMax,
        probe.protonMeanMin,
        probe.protonMeanMax,
        true,
        timingSigmaMin,
        probe.sigmaMax,
        probe.sigmaInitial,
        minimumGlobalSeparation,
        minimumGlobalAmplitudeSignificance,
        useDeviancePerEntryValidation,
        maximumGlobalPoissonDevianceNdf,
        maximumGlobalPoissonDeviancePerEntry,
        globalBoundFractionTolerance,
        minimumGlobalSliceEntries
      );

      if (shape.valid) {
        ++validShapes;
        validSeparationSum += shape.separation;
        validDevianceNdfSum += shape.poissonDevianceNdf;
      }

      if (shape.fitFunction) {
        delete shape.fitFunction;
        shape.fitFunction = nullptr;
      }

      delete projection;
    }

    delete probePID;
    probe.validShapes = validShapes;

    if (validShapes > 0) {
      probe.meanSeparation =
        validSeparationSum / validShapes;
      probe.meanPoissonDevianceNdf =
        validDevianceNdfSum / validShapes;
    }

    return validShapes;
  };

  auto countValidGlobalShapesForCT = [
    &
  ](
    TimingBranchProbe &probe
  ) -> int {
    probe.branch = ctTimeBranch;
    probe.displayMin = ctTimeMin;
    probe.displayMax = ctTimeMax;
    probe.fitMin = ctTimingFitMin;
    probe.fitMax = ctTimingFitMax;
    probe.kaonMeanMin = ctKaonMeanMin;
    probe.kaonMeanMax = ctKaonMeanMax;
    probe.protonMeanMin = ctProtonMeanMin;
    probe.protonMeanMax = ctProtonMeanMax;
    probe.sigmaMax = timingSigmaMax;
    probe.sigmaInitial = timingSigmaInitial;
    probe.histogramBins = ctTimeHistogramBins;

    const std::string ctRangeCut =
      TString::Format(
        "%s >= %.17g && %s <= %.17g && "
        "%s >= %.17g && %s <= %.17g",
        aeroBranch.c_str(),
        aeroMin,
        aeroBranch.c_str(),
        aeroMax,
        ctTimeBranch.c_str(),
        probe.displayMin,
        ctTimeBranch.c_str(),
        probe.displayMax
      ).Data();

    const std::string ctAnalysisCut =
      "(" + acceptanceCut + ") && (" +
      ctRangeCut + ")";

    gROOT->cd();

    auto *probePID = new TH2D(
      "h_ct_timing_probe_pid",
      "CT global timing probe",
      nAeroHistogramBins,
      aeroMin,
      aeroMax,
      probe.histogramBins,
      probe.displayMin,
      probe.displayMax
    );

    probePID->Sumw2();

    tree->Draw(
      TString::Format(
        "%s:%s>>%s",
        ctTimeBranch.c_str(),
        aeroBranch.c_str(),
        probePID->GetName()
      ),
      ctAnalysisCut.c_str(),
      "goff"
    );

    probePID->SetDirectory(nullptr);

    int validShapes = 0;
    double validSeparationSum = 0.0;
    double validDevianceNdfSum = 0.0;

    for (
      int aeroSlice = 0;
      aeroSlice < nAeroSlices;
      ++aeroSlice
    ) {
      const double aeroLow =
        aeroEdges.at(aeroSlice);

      const double aeroHigh =
        aeroEdges.at(aeroSlice + 1);

      const int firstXBin = std::max(
        1,
        probePID
          ->GetXaxis()
          ->FindFixBin(
            std::nextafter(
              aeroLow,
              aeroHigh
            )
          )
      );

      const int lastXBin = std::min(
        probePID->GetNbinsX(),
        probePID
          ->GetXaxis()
          ->FindFixBin(
            std::nextafter(
              aeroHigh,
              aeroLow
            )
          )
      );

      auto *projection = probePID->ProjectionY(
        TString::Format(
          "h_ct_timing_probe_aero_%d",
          aeroSlice
        ),
        firstXBin,
        lastXBin,
        "e"
      );

      projection->SetDirectory(nullptr);

      TimingShape shape = fitGlobalTimingShape(
        projection,
        TString::Format(
          "f_ct_timing_probe_aero_%d",
          aeroSlice
        ).Data(),
        probe.fitMin,
        probe.fitMax,
        probe.kaonMeanMin,
        probe.kaonMeanMax,
        probe.protonMeanMin,
        probe.protonMeanMax,
        false,
        timingSigmaMin,
        probe.sigmaMax,
        probe.sigmaInitial,
        minimumGlobalSeparation,
        minimumGlobalAmplitudeSignificance,
        useDeviancePerEntryValidation,
        maximumGlobalPoissonDevianceNdf,
        maximumGlobalPoissonDeviancePerEntry,
        globalBoundFractionTolerance,
        minimumGlobalSliceEntries
      );

      if (shape.valid) {
        ++validShapes;
        validSeparationSum += shape.separation;
        validDevianceNdfSum += shape.poissonDevianceNdf;
      }

      if (shape.fitFunction) {
        delete shape.fitFunction;
        shape.fitFunction = nullptr;
      }

      delete projection;
    }

    delete probePID;
    probe.validShapes = validShapes;

    if (validShapes > 0) {
      probe.meanSeparation =
        validSeparationSum / validShapes;
      probe.meanPoissonDevianceNdf =
        validDevianceNdfSum / validShapes;
    }

    return validShapes;
  };

  std::vector<TimingBranchProbe> rfProbes;

  for (const std::string &candidate : rfBranchCandidates) {
    if (
      candidate.empty() ||
      !tree->GetBranch(candidate.c_str())
    ) {
      continue;
    }

    rfTimingAttempted = true;

    TimingBranchProbe probe;

    std::cout
      << "RF timing probe using branch "
      << candidate
      << " from tree "
      << treeName
      << std::endl;

    const int validShapes =
      countValidGlobalShapesForBranch(
        candidate,
        probe
      );

    std::cout
      << "  valid global shapes: "
      << validShapes
      << " / "
      << nAeroSlices
      << "; range ["
      << probe.displayMin
      << ", "
      << probe.displayMax
      << "]; fit ["
      << probe.fitMin
      << ", "
      << probe.fitMax
      << "]; bins "
      << probe.histogramBins
      << "; mean separation "
      << probe.meanSeparation
      << "; mean D/ndf "
      << probe.meanPoissonDevianceNdf
      << "; peak seeds "
      << (
        probe.peakPairFound
          ? TString::Format(
              "p=%.4g, K=%.4g",
              probe.protonSeedMean,
              probe.kaonSeedMean
            ).Data()
          : "not resolved"
      )
      << std::endl;

    rfProbes.push_back(probe);
  }

  TimingBranchProbe ctProbe;
  ctTimingEvaluated = true;

  std::cout
    << "CT timing probe using branch "
    << ctTimeBranch
    << " from tree "
    << treeName
    << std::endl;

  ctProbeValidShapes =
    countValidGlobalShapesForCT(ctProbe);

  std::cout
    << "  valid global shapes: "
    << ctProbe.validShapes
    << " / "
    << nAeroSlices
    << "; range ["
    << ctProbe.displayMin
    << ", "
    << ctProbe.displayMax
    << "]; fit ["
    << ctProbe.fitMin
    << ", "
    << ctProbe.fitMax
    << "]; bins "
    << ctProbe.histogramBins
    << "; mean separation "
    << ctProbe.meanSeparation
    << "; mean D/ndf "
    << ctProbe.meanPoissonDevianceNdf
    << std::endl;

  const TimingBranchProbe *bestRFProbe = nullptr;

  if (!rfProbes.empty()) {
    const auto bestProbeIter = std::max_element(
      rfProbes.begin(),
      rfProbes.end(),
      [](
        const TimingBranchProbe &left,
        const TimingBranchProbe &right
      ) {
        if (left.validShapes != right.validShapes) {
          return left.validShapes < right.validShapes;
        }

        if (
          std::abs(
            left.meanSeparation -
            right.meanSeparation
          ) > 1.0e-9
        ) {
          return left.meanSeparation < right.meanSeparation;
        }

        return
          left.meanPoissonDevianceNdf >
          right.meanPoissonDevianceNdf;
      }
    );

    if (bestProbeIter != rfProbes.end()) {
      bestRFProbe = &(*bestProbeIter);
      rfProbeValidShapes = bestRFProbe->validShapes;
    }
  }

  auto compareTimingProbes = [](
    const TimingBranchProbe &left,
    const TimingBranchProbe &right
  ) -> int {
    if (left.validShapes != right.validShapes) {
      return left.validShapes > right.validShapes ? 1 : -1;
    }

    if (
      std::abs(
        left.meanSeparation -
        right.meanSeparation
      ) > 1.0e-9
    ) {
      return left.meanSeparation > right.meanSeparation ? 1 : -1;
    }

    const bool leftFinite =
      std::isfinite(left.meanPoissonDevianceNdf);

    const bool rightFinite =
      std::isfinite(right.meanPoissonDevianceNdf);

    if (leftFinite != rightFinite) {
      return leftFinite ? 1 : -1;
    }

    if (
      leftFinite &&
      std::abs(
        left.meanPoissonDevianceNdf -
        right.meanPoissonDevianceNdf
      ) > 1.0e-9
    ) {
      return
        left.meanPoissonDevianceNdf <
        right.meanPoissonDevianceNdf
          ? 1
          : -1;
    }

    return 0;
  };

  const bool selectRF =
    bestRFProbe &&
    bestRFProbe->validShapes > 0 &&
    compareTimingProbes(*bestRFProbe, ctProbe) >= 0;

  if (selectRF) {
    rfTimeBranch = bestRFProbe->branch;
    timeBranch = rfTimeBranch;
    rfTimingSelected = true;

    timeMin = bestRFProbe->displayMin;
    timeMax = bestRFProbe->displayMax;
    timingFitMin = bestRFProbe->fitMin;
    timingFitMax = bestRFProbe->fitMax;

    kaonMeanMin = bestRFProbe->kaonMeanMin;
    kaonMeanMax = bestRFProbe->kaonMeanMax;
    protonMeanMin = bestRFProbe->protonMeanMin;
    protonMeanMax = bestRFProbe->protonMeanMax;
    timingSigmaMax = bestRFProbe->sigmaMax;
    timingSigmaInitial = bestRFProbe->sigmaInitial;
    nTimeHistogramBins = bestRFProbe->histogramBins;

    timingSelectionReason =
      compareTimingProbes(*bestRFProbe, ctProbe) > 0
        ? "rf_probe_ranked_better_than_ct"
        : "rf_won_exact_probe_tie";

    std::cout
      << "Selected RF timing branch: "
      << rfTimeBranch
      << " (RF valid="
      << bestRFProbe->validShapes
      << ", mean separation="
      << bestRFProbe->meanSeparation
      << ", mean D/ndf="
      << bestRFProbe->meanPoissonDevianceNdf
      << "; CT valid="
      << ctProbe.validShapes
      << ", mean separation="
      << ctProbe.meanSeparation
      << ", mean D/ndf="
      << ctProbe.meanPoissonDevianceNdf
      << ")"
      << std::endl;
  } else {
    timeBranch = ctTimeBranch;
    rfTimingSelected = false;
    ctFallbackUsed = false;

    timeMin = ctProbe.displayMin;
    timeMax = ctProbe.displayMax;
    timingFitMin = ctProbe.fitMin;
    timingFitMax = ctProbe.fitMax;

    kaonMeanMin = ctProbe.kaonMeanMin;
    kaonMeanMax = ctProbe.kaonMeanMax;
    protonMeanMin = ctProbe.protonMeanMin;
    protonMeanMax = ctProbe.protonMeanMax;
    timingSigmaMax = ctProbe.sigmaMax;
    timingSigmaInitial = ctProbe.sigmaInitial;
    nTimeHistogramBins = ctProbe.histogramBins;

    if (!bestRFProbe) {
      timingSelectionReason = disableRFTiming
        ? "rf_disabled_ct_selected"
        : "no_rf_branch_ct_selected";
    } else if (ctProbe.validShapes == 0 && bestRFProbe->validShapes == 0) {
      timingSelectionReason =
        "no_valid_rf_or_ct_shapes_ct_selected_for_diagnostics";
    } else {
      timingSelectionReason =
        "ct_probe_ranked_better_than_rf";
    }

    std::cout
      << "Selected CT timing branch: "
      << ctTimeBranch
      << " (CT valid="
      << ctProbe.validShapes
      << ", mean separation="
      << ctProbe.meanSeparation
      << ", mean D/ndf="
      << ctProbe.meanPoissonDevianceNdf;

    if (bestRFProbe) {
      std::cout
        << "; best RF "
        << bestRFProbe->branch
        << " valid="
        << bestRFProbe->validShapes
        << ", mean separation="
        << bestRFProbe->meanSeparation
        << ", mean D/ndf="
        << bestRFProbe->meanPoissonDevianceNdf;
    }

    std::cout
      << ")"
      << std::endl;

    if (!bestRFProbe && !disableRFTiming) {
      printRFLikeBranches(tree);
    }
  }

  // RF: proton is the left/lower-time peak.
  // CTime_ROC1: proton is the right/higher-time peak.
  const bool protonPeakIsLower = rfTimingSelected;

  const std::string timingAxisTitle =
    timeBranch + " [ns]";

  std::cout
    << "Selected timing branch: "
    << timeBranch
    << " for "
    << epsSetting
    << " epsilon; display range ["
    << timeMin
    << ", "
    << timeMax
    << "] ns with "
    << nTimeHistogramBins
    << " bins; beam-bunch spacing = "
    << beamBunchSpacingNs
    << " ns; proton peak expected on the "
    << (
      protonPeakIsLower
        ? "left/lower-time side"
        : "right/higher-time side"
    )
    << "; RF phase correction/wrapping = disabled"
    << std::endl;

  const std::string pidRangeCut =
    TString::Format(
      "%s >= %.17g && %s <= %.17g && "
      "%s >= %.17g && %s <= %.17g",
      aeroBranch.c_str(),
      aeroMin,
      aeroBranch.c_str(),
      aeroMax,
      timeBranch.c_str(),
      timeMin,
      timeBranch.c_str(),
      timeMax
    ).Data();

  const std::string analysisCut =
    "(" + acceptanceCut + ") && (" +
    pidRangeCut + ")";

  // ------------------------------------------------------------------
  // Global PID histogram.
  // ------------------------------------------------------------------

  gROOT->cd();

  auto *hGlobalPID = new TH2D(
    "h_global_pid",
    TString::Format(
      "Global %s vs P_aero_npeSum;"
      "P_aero_npeSum;%s;Counts",
      timeBranch.c_str(),
      timingAxisTitle.c_str()
    ),
    nAeroHistogramBins,
    aeroMin,
    aeroMax,
    nTimeHistogramBins,
    timeMin,
    timeMax
  );

  hGlobalPID->Sumw2();

  tree->Draw(
    TString::Format(
      "%s:%s>>%s",
      timeBranch.c_str(),
      aeroBranch.c_str(),
      hGlobalPID->GetName()
    ),
    analysisCut.c_str(),
    "goff"
  );

  hGlobalPID->SetDirectory(nullptr);

  if (hGlobalPID->Integral() <= 0.0) {
    std::cerr
      << "Global PID histogram is empty."
      << std::endl;

    inputFile->Close();
    return;
  }

  std::cout
    << "Timing fit range: ["
    << timingFitMin
    << ", "
    << timingFitMax
    << "] ns using "
    << timeBranch
    << std::endl;

  // ------------------------------------------------------------------
  // Global timing fits.
  // ------------------------------------------------------------------

  std::vector<TH1D *> globalTimingProjections(
    nAeroSlices,
    nullptr
  );

  std::vector<TimingShape> globalShapes(
    nAeroSlices
  );

  for (
    int aeroSlice = 0;
    aeroSlice < nAeroSlices;
    ++aeroSlice
  ) {
    const double aeroLow =
      aeroEdges.at(aeroSlice);

    const double aeroHigh =
      aeroEdges.at(aeroSlice + 1);

    const int firstXBin = std::max(
      1,
      hGlobalPID
        ->GetXaxis()
        ->FindFixBin(
          std::nextafter(
            aeroLow,
            aeroHigh
          )
        )
    );

    const int lastXBin = std::min(
      hGlobalPID->GetNbinsX(),
      hGlobalPID
        ->GetXaxis()
        ->FindFixBin(
          std::nextafter(
            aeroHigh,
            aeroLow
          )
        )
    );

    auto *projection =
      hGlobalPID->ProjectionY(
        TString::Format(
          "h_global_time_aero_slice_%d",
          aeroSlice
        ),
        firstXBin,
        lastXBin,
        "e"
      );

    projection->SetDirectory(nullptr);

    projection->SetTitle(
      TString::Format(
        "Global timing fit: %.1f #leq aero < %.1f;"
        "%s;Counts",
        aeroLow,
        aeroHigh,
        timingAxisTitle.c_str()
      )
    );

    globalTimingProjections.at(aeroSlice) =
      projection;

    globalShapes.at(aeroSlice) =
      fitGlobalTimingShape(
        projection,
        TString::Format(
          "f_global_time_aero_slice_%d",
          aeroSlice
        ).Data(),
        timingFitMin,
        timingFitMax,
        kaonMeanMin,
        kaonMeanMax,
        protonMeanMin,
        protonMeanMax,
        protonPeakIsLower,
        timingSigmaMin,
        timingSigmaMax,
        timingSigmaInitial,
        minimumGlobalSeparation,
        minimumGlobalAmplitudeSignificance,
        useDeviancePerEntryValidation,
        maximumGlobalPoissonDevianceNdf,
        maximumGlobalPoissonDeviancePerEntry,
        globalBoundFractionTolerance,
        minimumGlobalSliceEntries
      );
  }

  int validGlobalShapes = 0;

  for (
    const TimingShape &shape :
    globalShapes
  ) {
    if (shape.valid) {
      ++validGlobalShapes;
    }
  }

  if (validGlobalShapes == 0) {
    std::cerr
      << "No identifiable proton-kaon timing shapes were found. "
      << "[v9 diagnostic mode] Writing raw diagnostic plots to "
      << outputPDF
      << " and "
      << outputROOT
      << "."
      << std::endl;

    // --------------------------------------------------------------
    // Diagnostic-only output.  Do not attempt proton subtraction when
    // none of the global timing fits passes the shape-quality cuts.
    // Still save the raw global and per-delta PID distributions so the
    // failed fits and underlying timing structure can be inspected.
    // --------------------------------------------------------------

    auto *diagnosticGlobalCanvas = new TCanvas(
      "canvas_global_timing_diagnostics",
      "Global timing diagnostics",
      2100,
      1300
    );

    diagnosticGlobalCanvas->Divide(3, 2);
    diagnosticGlobalCanvas->cd(1);
    gPad->SetRightMargin(0.16);
    gPad->SetLogz();
    hGlobalPID->Draw("COLZ");

    auto *globalWarning = new TPaveText(
      0.12,
      0.78,
      0.82,
      0.90,
      "NDC"
    );

    globalWarning->SetFillColor(kWhite);
    globalWarning->SetFillStyle(1001);
    globalWarning->SetBorderSize(1);
    globalWarning->SetTextColor(kRed + 1);
    globalWarning->SetTextAlign(12);
    globalWarning->AddText(
      "Diagnostic only: no global timing shape passed validation"
    );
    globalWarning->Draw();

    for (
      int aeroSlice = 0;
      aeroSlice < nAeroSlices;
      ++aeroSlice
    ) {
      diagnosticGlobalCanvas->cd(aeroSlice + 2);

      TH1D *projection =
        globalTimingProjections.at(aeroSlice);

      if (!projection) {
        continue;
      }

      projection->SetLineColor(kBlack);
      projection->SetLineWidth(2);
      projection->Draw("E");

      const TimingShape &shape =
        globalShapes.at(aeroSlice);

      if (shape.fitFunction) {
        shape.fitFunction->SetLineColor(kOrange + 7);
        shape.fitFunction->SetLineWidth(2);
        shape.fitFunction->Draw("SAME");
      }

      auto *fitStatusText = new TPaveText(
        0.49,
        0.51,
        0.89,
        0.89,
        "NDC"
      );

      fitStatusText->SetFillStyle(0);
      fitStatusText->SetBorderSize(0);
      fitStatusText->SetTextAlign(12);
      fitStatusText->AddText(
        TString::Format(
          "entries: %.0f",
          projection->GetEntries()
        )
      );
      fitStatusText->AddText(
        TString::Format(
          "fit status: %d",
          shape.fitStatus
        )
      );
      fitStatusText->AddText(
        TString::Format(
          "bound hit: %s",
          shape.boundHit ? "yes" : "no"
        )
      );
      fitStatusText->AddText(
        TString::Format(
          "K #mu/#sigma: %.3f / %.3f",
          shape.kaonMean,
          shape.kaonSigma
        )
      );
      fitStatusText->AddText(
        TString::Format(
          "p #mu/#sigma: %.3f / %.3f",
          shape.protonMean,
          shape.protonSigma
        )
      );
      fitStatusText->AddText(
        TString::Format(
          "separation: %.2f",
          shape.separation
        )
      );
      fitStatusText->AddText(
        TString::Format(
          "K/p significance: %.1f / %.1f",
          shape.kaonSignificance,
          shape.protonSignificance
        )
      );
      fitStatusText->AddText(
        TString::Format(
          "Poisson D/ndf: %.2f (%d)",
          shape.poissonDevianceNdf,
          shape.goodnessNdf
        )
      );
      fitStatusText->AddText(
        TString::Format(
          "D/N: %.3f",
          shape.poissonDeviancePerEntry
        )
      );
      fitStatusText->Draw();
    }

    std::vector<TH2D *> diagnosticDeltaPIDHistograms(
      nDeltaBins,
      nullptr
    );

    std::vector<std::vector<TH1D *>>
      diagnosticDeltaTimingProjections(
        nDeltaBins,
        std::vector<TH1D *>(
          nAeroSlices,
          nullptr
        )
      );

    std::vector<TCanvas *> diagnosticDeltaCanvases(
      nDeltaBins,
      nullptr
    );

    for (
      int deltaBin = 0;
      deltaBin < nDeltaBins;
      ++deltaBin
    ) {
      const double deltaLow =
        deltaMin + deltaBin * deltaBinWidth;

      const double deltaHigh =
        deltaLow + deltaBinWidth;

      const std::string deltaCut =
        TString::Format(
          "%s >= %.17g && %s %s %.17g",
          deltaBranch.c_str(),
          deltaLow,
          deltaBranch.c_str(),
          deltaBin == nDeltaBins - 1 ? "<=" : "<",
          deltaHigh
        ).Data();

      const std::string fullDeltaCut =
        "(" + analysisCut + ") && (" +
        deltaCut + ")";

      gROOT->cd();

      auto *deltaPID = new TH2D(
        TString::Format(
          "h_diagnostic_pid_delta_bin_%d",
          deltaBin
        ),
        TString::Format(
          "Raw PID diagnostic: %.1f #leq #delta %s %.1f;"
          "P_aero_npeSum;%s;Counts",
          deltaLow,
          deltaBin == nDeltaBins - 1 ? "#leq" : "<",
          deltaHigh,
          timingAxisTitle.c_str()
        ),
        nAeroHistogramBins,
        aeroMin,
        aeroMax,
        nTimeHistogramBins,
        timeMin,
        timeMax
      );

      deltaPID->Sumw2();

      tree->Draw(
        TString::Format(
          "%s:%s>>%s",
          timeBranch.c_str(),
          aeroBranch.c_str(),
          deltaPID->GetName()
        ),
        fullDeltaCut.c_str(),
        "goff"
      );

      deltaPID->SetDirectory(nullptr);
      diagnosticDeltaPIDHistograms.at(deltaBin) =
        deltaPID;

      auto *deltaCanvas = new TCanvas(
        TString::Format(
          "canvas_raw_delta_diagnostic_%d",
          deltaBin
        ),
        TString::Format(
          "Raw delta diagnostic %.1f to %.1f",
          deltaLow,
          deltaHigh
        ),
        2100,
        1300
      );

      deltaCanvas->Divide(3, 2);
      deltaCanvas->cd(1);
      gPad->SetRightMargin(0.16);

      if (deltaPID->Integral() > 0.0) {
        gPad->SetLogz();
      }

      deltaPID->Draw("COLZ");

      auto *deltaLabel = new TPaveText(
        0.12,
        0.78,
        0.75,
        0.90,
        "NDC"
      );

      deltaLabel->SetFillColor(kWhite);
      deltaLabel->SetFillStyle(1001);
      deltaLabel->SetBorderSize(1);
      deltaLabel->SetTextAlign(12);
      deltaLabel->AddText(
        TString::Format(
          "raw entries: %.0f",
          deltaPID->GetEntries()
        )
      );
      deltaLabel->AddText(
        "No proton-cleaning weights applied"
      );
      deltaLabel->Draw();

      for (
        int aeroSlice = 0;
        aeroSlice < nAeroSlices;
        ++aeroSlice
      ) {
        const double aeroLow =
          aeroEdges.at(aeroSlice);

        const double aeroHigh =
          aeroEdges.at(aeroSlice + 1);

        const int firstXBin = std::max(
          1,
          deltaPID
            ->GetXaxis()
            ->FindFixBin(
              std::nextafter(
                aeroLow,
                aeroHigh
              )
            )
        );

        const int lastXBin = std::min(
          deltaPID->GetNbinsX(),
          deltaPID
            ->GetXaxis()
            ->FindFixBin(
              std::nextafter(
                aeroHigh,
                aeroLow
              )
            )
        );

        auto *projection =
          deltaPID->ProjectionY(
            TString::Format(
              "h_diagnostic_time_delta_%d_aero_%d",
              deltaBin,
              aeroSlice
            ),
            firstXBin,
            lastXBin,
            "e"
          );

        projection->SetDirectory(nullptr);
        projection->SetTitle(
          TString::Format(
            "%.1f #leq #delta %s %.1f, "
            "%.1f #leq aero < %.1f;"
            "%s;Counts",
            deltaLow,
            deltaBin == nDeltaBins - 1 ? "#leq" : "<",
            deltaHigh,
            aeroLow,
            aeroHigh,
            timingAxisTitle.c_str()
          )
        );
        projection->SetLineColor(kBlack);
        projection->SetLineWidth(2);

        diagnosticDeltaTimingProjections
          .at(deltaBin)
          .at(aeroSlice) = projection;

        deltaCanvas->cd(aeroSlice + 2);
        projection->Draw("E");

        auto *entryLabel = new TPaveText(
          0.58,
          0.78,
          0.88,
          0.88,
          "NDC"
        );

        entryLabel->SetFillStyle(0);
        entryLabel->SetBorderSize(0);
        entryLabel->SetTextAlign(12);
        entryLabel->AddText(
          TString::Format(
            "entries: %.0f",
            projection->GetEntries()
          )
        );
        entryLabel->Draw();
      }

      diagnosticDeltaCanvases.at(deltaBin) =
        deltaCanvas;
    }

    diagnosticGlobalCanvas->Modified();
    diagnosticGlobalCanvas->Update();
    diagnosticGlobalCanvas->Print(
      (outputPDF + "[").c_str()
    );
    diagnosticGlobalCanvas->Print(
      outputPDF.c_str()
    );

    for (
      TCanvas *canvas :
      diagnosticDeltaCanvases
    ) {
      if (!canvas) {
        continue;
      }

      canvas->Modified();
      canvas->Update();
      canvas->Print(outputPDF.c_str());
    }

    diagnosticGlobalCanvas->Print(
      (outputPDF + "]").c_str()
    );

    TFile diagnosticOutputFile(
      outputROOT.c_str(),
      "RECREATE"
    );

    if (diagnosticOutputFile.IsZombie()) {
      std::cerr
        << "Unable to create diagnostic output file: "
        << outputROOT
        << std::endl;

      inputFile->Close();
      return;
    }

    TDirectory *diagnosticGlobalDirectory =
      diagnosticOutputFile.mkdir("global_diagnostics");

    diagnosticGlobalDirectory->cd();
    hGlobalPID->Write();
    diagnosticGlobalCanvas->Write();

    for (
      int aeroSlice = 0;
      aeroSlice < nAeroSlices;
      ++aeroSlice
    ) {
      if (globalTimingProjections.at(aeroSlice)) {
        globalTimingProjections
          .at(aeroSlice)
          ->Write();
      }

      if (globalShapes.at(aeroSlice).fitFunction) {
        globalShapes
          .at(aeroSlice)
          .fitFunction
          ->Write();
      }
    }

    diagnosticOutputFile.cd();

    TDirectory *diagnosticDeltaDirectory =
      diagnosticOutputFile.mkdir("delta_diagnostics");

    for (
      int deltaBin = 0;
      deltaBin < nDeltaBins;
      ++deltaBin
    ) {
      diagnosticDeltaDirectory->cd();

      TDirectory *deltaDirectory =
        diagnosticDeltaDirectory->mkdir(
          TString::Format(
            "delta_bin_%02d",
            deltaBin
          )
        );

      deltaDirectory->cd();

      if (diagnosticDeltaPIDHistograms.at(deltaBin)) {
        diagnosticDeltaPIDHistograms
          .at(deltaBin)
          ->Write();
      }

      for (
        int aeroSlice = 0;
        aeroSlice < nAeroSlices;
        ++aeroSlice
      ) {
        if (
          diagnosticDeltaTimingProjections
            .at(deltaBin)
            .at(aeroSlice)
        ) {
          diagnosticDeltaTimingProjections
            .at(deltaBin)
            .at(aeroSlice)
            ->Write();
        }
      }

      if (diagnosticDeltaCanvases.at(deltaBin)) {
        diagnosticDeltaCanvases
          .at(deltaBin)
          ->Write();
      }
    }

    diagnosticOutputFile.cd();

    TNamed diagnosticMacroVersion(
      "macro_version",
      macroVersion
    );
    diagnosticMacroVersion.Write();

    TNamed diagnosticMode(
      "analysis_mode",
      "diagnostic_only_no_valid_global_timing_shapes"
    );
    diagnosticMode.Write();

    TParameter<int> numberOfValidGlobalShapes(
      "valid_global_timing_shapes",
      validGlobalShapes
    );
    numberOfValidGlobalShapes.Write();

    TNamed diagnosticTimingBranch(
      "timing_branch",
      timeBranch.c_str()
    );
    diagnosticTimingBranch.Write();

    TNamed diagnosticRequestedRFBranch(
      "rf_timing_branch",
      rfTimeBranch.c_str()
    );
    diagnosticRequestedRFBranch.Write();

    TParameter<int>(
      "rf_timing_attempted",
      rfTimingAttempted ? 1 : 0
    ).Write();

    TParameter<int>(
      "rf_timing_selected",
      rfTimingSelected ? 1 : 0
    ).Write();

    TParameter<int>(
      "proton_peak_is_lower",
      protonPeakIsLower ? 1 : 0
    ).Write();

    TParameter<int>(
      "ct_fallback_used",
      ctFallbackUsed ? 1 : 0
    ).Write();

    TParameter<int>(
      "ct_timing_evaluated",
      ctTimingEvaluated ? 1 : 0
    ).Write();

    TParameter<int>(
      "rf_probe_valid_global_shapes",
      rfProbeValidShapes
    ).Write();

    TParameter<int>(
      "ct_probe_valid_global_shapes",
      ctProbeValidShapes
    ).Write();

    TNamed diagnosticSelectionReason(
      "timing_selection_reason",
      timingSelectionReason.c_str()
    );
    diagnosticSelectionReason.Write();

    TNamed diagnosticFitStatistic(
      "fit_statistic",
      "Baker-Cousins Poisson deviance; bin-integrated likelihood"
    );
    diagnosticFitStatistic.Write();

    TNamed diagnosticValidationStatistic(
      "fit_validation_statistic",
      useDeviancePerEntryValidation
        ? "poisson_deviance_per_entry"
        : "poisson_deviance_per_ndf"
    );
    diagnosticValidationStatistic.Write();

    TParameter<double>(
      "maximum_global_poisson_deviance_per_ndf",
      maximumGlobalPoissonDevianceNdf
    ).Write();

    TParameter<double>(
      "maximum_global_poisson_deviance_per_entry",
      maximumGlobalPoissonDeviancePerEntry
    ).Write();

    TParameter<int>(
      "rf_correction_used",
      0
    ).Write();

    TParameter<double>(
      "beam_bunch_spacing_ns",
      beamBunchSpacingNs
    ).Write();

    TParameter<double>(
      "timing_fit_min_ns",
      timingFitMin
    ).Write();

    TParameter<double>(
      "timing_fit_max_ns",
      timingFitMax
    ).Write();

    diagnosticOutputFile.Write();
    diagnosticOutputFile.Close();

    std::cout
      << "Diagnostic PDF written to: "
      << outputPDF
      << std::endl;

    std::cout
      << "Diagnostic ROOT file written to: "
      << outputROOT
      << std::endl;

    inputFile->Close();
    return;
  }

  // ------------------------------------------------------------------
  // Per-delta and per-aerogel fits.
  // ------------------------------------------------------------------

  std::vector<TH2D *> deltaPIDHistograms(
    nDeltaBins,
    nullptr
  );

  std::vector<std::vector<TH1D *>>
    deltaTimingProjections(
      nDeltaBins,
      std::vector<TH1D *>(
        nAeroSlices,
        nullptr
      )
    );

  std::vector<std::vector<SliceFitResult>>
    deltaSliceFits(
      nDeltaBins,
      std::vector<SliceFitResult>(
        nAeroSlices
      )
    );

  std::vector<double> protonYieldByDelta(
    nDeltaBins,
    0.0
  );

  std::vector<double> kaonYieldByDelta(
    nDeltaBins,
    0.0
  );

  std::vector<double> otherYieldByDelta(
    nDeltaBins,
    0.0
  );

  std::vector<double> dataYieldByDelta(
    nDeltaBins,
    0.0
  );

  std::vector<double> fittedDataYieldByDelta(
    nDeltaBins,
    0.0
  );

  std::vector<double> modelYieldByDelta(
    nDeltaBins,
    0.0
  );

  std::vector<double> poissonDevianceNdfByDelta(
    nDeltaBins,
    0.0
  );

  std::vector<int> validSlicesByDelta(
    nDeltaBins,
    0
  );

  std::vector<double> validCoverageByDelta(
    nDeltaBins,
    0.0
  );

  std::vector<SupportClass> supportByDelta(
    nDeltaBins,
    SupportClass::Unsupported
  );

  for (
    int deltaBin = 0;
    deltaBin < nDeltaBins;
    ++deltaBin
  ) {
    const double deltaLow =
      deltaMin +
      deltaBin * deltaBinWidth;

    const double deltaHigh =
      deltaLow + deltaBinWidth;

    const bool includeUpperEdge =
      deltaBin == nDeltaBins - 1;

    const std::string deltaCut =
      makeDeltaCut(
        deltaBranch,
        deltaLow,
        deltaHigh,
        includeUpperEdge
      );

    const std::string fullCut =
      "(" + analysisCut + ") && (" +
      deltaCut + ")";

    gROOT->cd();

    auto *hPID = new TH2D(
      TString::Format(
        "h_pid_delta_%d",
        deltaBin
      ),
      TString::Format(
        "PID plane: %.1f #leq #delta %s %.1f;"
        "P_aero_npeSum;%s;Counts",
        deltaLow,
        includeUpperEdge ? "#leq" : "<",
        deltaHigh,
        timingAxisTitle.c_str()
      ),
      nAeroHistogramBins,
      aeroMin,
      aeroMax,
      nTimeHistogramBins,
      timeMin,
      timeMax
    );

    hPID->Sumw2();

    tree->Draw(
      TString::Format(
        "%s:%s>>%s",
        timeBranch.c_str(),
        aeroBranch.c_str(),
        hPID->GetName()
      ),
      fullCut.c_str(),
      "goff"
    );

    hPID->SetDirectory(nullptr);

    deltaPIDHistograms.at(deltaBin) =
      hPID;

    dataYieldByDelta.at(deltaBin) =
      hPID->Integral();

    double devianceWeightedSum = 0.0;
    double devianceWeight = 0.0;

    for (
      int aeroSlice = 0;
      aeroSlice < nAeroSlices;
      ++aeroSlice
    ) {
      const double aeroLow =
        aeroEdges.at(aeroSlice);

      const double aeroHigh =
        aeroEdges.at(aeroSlice + 1);

      const int firstXBin = std::max(
        1,
        hPID
          ->GetXaxis()
          ->FindFixBin(
            std::nextafter(
              aeroLow,
              aeroHigh
            )
          )
      );

      const int lastXBin = std::min(
        hPID->GetNbinsX(),
        hPID
          ->GetXaxis()
          ->FindFixBin(
            std::nextafter(
              aeroHigh,
              aeroLow
            )
          )
      );

      auto *projection =
        hPID->ProjectionY(
          TString::Format(
            "h_time_delta_%d_aero_%d",
            deltaBin,
            aeroSlice
          ),
          firstXBin,
          lastXBin,
          "e"
        );

      projection->SetDirectory(nullptr);

      projection->SetTitle(
        TString::Format(
          "%.1f #leq #delta %s %.1f, "
          "%.1f #leq aero < %.1f;"
          "%s;Counts",
          deltaLow,
          includeUpperEdge ? "#leq" : "<",
          deltaHigh,
          aeroLow,
          aeroHigh,
          timingAxisTitle.c_str()
        )
      );

      deltaTimingProjections
        .at(deltaBin)
        .at(aeroSlice) =
          projection;

      if (!globalShapes.at(aeroSlice).valid) {
        continue;
      }

      const SliceFitResult sliceFit =
        fitDeltaTimingSlice(
          projection,
          globalShapes.at(aeroSlice),
          TString::Format(
            "f_time_delta_%d_aero_%d",
            deltaBin,
            aeroSlice
          ).Data(),
          timingFitMin,
          timingFitMax,
          useDeviancePerEntryValidation,
          maximumSlicePoissonDevianceNdf,
          maximumSlicePoissonDeviancePerEntry,
          minimumSliceModelDataRatio,
          maximumSliceModelDataRatio,
          minimumDeltaSliceEntries
        );

      deltaSliceFits
        .at(deltaBin)
        .at(aeroSlice) =
          sliceFit;

      if (!sliceFit.valid) {
        continue;
      }

      protonYieldByDelta.at(deltaBin) +=
        sliceFit.protonYield;

      kaonYieldByDelta.at(deltaBin) +=
        sliceFit.kaonYield;

      otherYieldByDelta.at(deltaBin) +=
        sliceFit.otherYield;

      fittedDataYieldByDelta.at(deltaBin) +=
        sliceFit.dataYield;

      modelYieldByDelta.at(deltaBin) +=
        sliceFit.modelYield;

      ++validSlicesByDelta.at(deltaBin);

      devianceWeightedSum +=
        sliceFit.poissonDevianceNdf *
        sliceFit.dataYield;

      devianceWeight +=
        sliceFit.dataYield;
    }

    if (devianceWeight > 0.0) {
      poissonDevianceNdfByDelta.at(deltaBin) =
        devianceWeightedSum /
        devianceWeight;
    }

    if (dataYieldByDelta.at(deltaBin) > 0.0) {
      validCoverageByDelta.at(deltaBin) =
        fittedDataYieldByDelta.at(deltaBin) /
        dataYieldByDelta.at(deltaBin);
    }

    const double totalModeledYield =
      protonYieldByDelta.at(deltaBin) +
      kaonYieldByDelta.at(deltaBin) +
      otherYieldByDelta.at(deltaBin);

    if (
      validSlicesByDelta.at(deltaBin) >=
        minimumSupportedSlices &&
      validCoverageByDelta.at(deltaBin) >=
        minimumSupportedCoverage &&
      totalModeledYield >= minimumModeledYield
    ) {
      supportByDelta.at(deltaBin) =
        SupportClass::Supported;
    } else if (
      validSlicesByDelta.at(deltaBin) >=
        minimumMarginalSlices &&
      validCoverageByDelta.at(deltaBin) >=
        minimumMarginalCoverage &&
      totalModeledYield >= minimumModeledYield
    ) {
      supportByDelta.at(deltaBin) =
        SupportClass::Marginal;
    } else {
      supportByDelta.at(deltaBin) =
        SupportClass::Unsupported;
    }
  }

  // ------------------------------------------------------------------
  // Summary histograms.
  // ------------------------------------------------------------------

  auto *hProtonYield = new TH1D(
    "h_proton_yield_delta",
    "Fitted yields versus SHMS #delta;"
    "SHMS #delta [%];"
    "Fitted yield",
    nDeltaBins,
    deltaMin,
    deltaMax
  );

  auto *hKaonYield = new TH1D(
    "h_kaon_yield_delta",
    "Fitted yields versus SHMS #delta;"
    "SHMS #delta [%];"
    "Fitted yield",
    nDeltaBins,
    deltaMin,
    deltaMax
  );

  auto *hOtherYield = new TH1D(
    "h_other_yield_delta",
    "Fitted yields versus SHMS #delta;"
    "SHMS #delta [%];"
    "Fitted yield",
    nDeltaBins,
    deltaMin,
    deltaMax
  );

  auto *hFitProtonWeightDelta = new TH1D(
    "h_fit_proton_weight_delta",
    "Integrated fitted proton fraction;"
    "SHMS #delta [%];"
    "w_{p}^{fit}(#delta)",
    nDeltaBins,
    deltaMin,
    deltaMax
  );

  auto *hFitChi2 = new TH1D(
    "h_fit_chi2_delta",
    "Average timing-fit quality;"
    "SHMS #delta [%];"
    "Weighted average Poisson D/ndf",
    nDeltaBins,
    deltaMin,
    deltaMax
  );

  auto *hFitCoverage = new TH1D(
    "h_fit_coverage_delta",
    "Accepted PID-fit coverage;"
    "SHMS #delta [%];"
    "Accepted fit data / total data",
    nDeltaBins,
    deltaMin,
    deltaMax
  );

  auto *hValidSlices = new TH1D(
    "h_valid_slices_delta",
    "Valid aerogel slices per #delta bin;"
    "SHMS #delta [%];"
    "Valid slices",
    nDeltaBins,
    deltaMin,
    deltaMax
  );

  auto *hSupportClass = new TH1D(
    "h_support_class_delta",
    "Proton-weight support class;"
    "SHMS #delta [%];"
    "Support class",
    nDeltaBins,
    deltaMin,
    deltaMax
  );

  const std::vector<TH1D *> deltaSummaryHistograms = {
    hProtonYield,
    hKaonYield,
    hOtherYield,
    hFitProtonWeightDelta,
    hFitChi2,
    hFitCoverage,
    hValidSlices,
    hSupportClass
  };

  for (
    TH1D *histogram :
    deltaSummaryHistograms
  ) {
    histogram->SetDirectory(nullptr);
    histogram->Sumw2();
  }

  for (
    int deltaBin = 0;
    deltaBin < nDeltaBins;
    ++deltaBin
  ) {
    const int rootBin =
      deltaBin + 1;

    const double protonYield =
      protonYieldByDelta.at(deltaBin);

    const double kaonYield =
      kaonYieldByDelta.at(deltaBin);

    const double otherYield =
      otherYieldByDelta.at(deltaBin);

    const double totalYield =
      protonYield +
      kaonYield +
      otherYield;

    hProtonYield->SetBinContent(
      rootBin,
      protonYield
    );

    hKaonYield->SetBinContent(
      rootBin,
      kaonYield
    );

    hOtherYield->SetBinContent(
      rootBin,
      otherYield
    );

    hFitChi2->SetBinContent(
      rootBin,
      poissonDevianceNdfByDelta.at(deltaBin)
    );

    hFitCoverage->SetBinContent(
      rootBin,
      validCoverageByDelta.at(deltaBin)
    );

    hValidSlices->SetBinContent(
      rootBin,
      validSlicesByDelta.at(deltaBin)
    );

    hSupportClass->SetBinContent(
      rootBin,
      static_cast<int>(
        supportByDelta.at(deltaBin)
      )
    );

    double protonWeight = 0.0;
    double protonWeightError = 0.0;

    if (
      supportByDelta.at(deltaBin) !=
        SupportClass::Unsupported &&
      totalYield >= minimumModeledYield
    ) {
      protonWeight =
        protonYield /
        totalYield;

      protonWeightError =
        std::sqrt(
          std::max(
            protonWeight *
            (1.0 - protonWeight) /
            totalYield,
            0.0
          )
        );

      if (
        supportByDelta.at(deltaBin) ==
          SupportClass::Marginal
      ) {
        protonWeightError *=
          marginalUncertaintyInflation;
      }
    }

    hFitProtonWeightDelta->SetBinContent(
      rootBin,
      std::clamp(
        protonWeight,
        0.0,
        1.0
      )
    );

    hFitProtonWeightDelta->SetBinError(
      rootBin,
      protonWeightError
    );
  }

  // ------------------------------------------------------------------
  // Event-level histograms.
  // ------------------------------------------------------------------

  auto *hRawMM = new TH1D(
    "h_mm_raw",
    "Kaon MM before and after proton cleaning;"
    "MM [GeV];"
    "Weighted counts",
    nMMBins,
    mmMin,
    mmMax
  );

  auto *hProtonMM = new TH1D(
    "h_mm_estimated_proton",
    "Estimated proton contamination;"
    "MM [GeV];"
    "Weighted counts",
    nMMBins,
    mmMin,
    mmMax
  );

  auto *hCleanedMM = new TH1D(
    "h_mm_proton_cleaned",
    "Proton-cleaned kaon MM;"
    "MM [GeV];"
    "Weighted counts",
    nMMBins,
    mmMin,
    mmMax
  );

  auto *hAppliedWeightSumDelta = new TH1D(
    "h_applied_weight_sum_delta",
    "Applied proton-weight sum;"
    "SHMS #delta [%];"
    "Sum of weights",
    nDeltaBins,
    deltaMin,
    deltaMax
  );

  auto *hAppliedWeightCountDelta = new TH1D(
    "h_applied_weight_count_delta",
    "Applied proton-weight event count;"
    "SHMS #delta [%];"
    "Events",
    nDeltaBins,
    deltaMin,
    deltaMax
  );

  auto *hAppliedWeightDelta = new TH1D(
    "h_applied_weight_delta",
    "Mean applied event-level proton weight;"
    "SHMS #delta [%];"
    "#LTw_{p}^{event}#GT",
    nDeltaBins,
    deltaMin,
    deltaMax
  );

  auto *hAppliedWeightSumMap = new TH2D(
    "h_applied_weight_sum_map",
    "Applied proton-weight sum;"
    "SHMS #delta [%];"
    "P_aero_npeSum;"
    "Sum of weights",
    nDeltaBins,
    deltaMin,
    deltaMax,
    nAeroSlices,
    aeroEdges.data()
  );

  auto *hAppliedWeightCountMap = new TH2D(
    "h_applied_weight_count_map",
    "Applied proton-weight event count;"
    "SHMS #delta [%];"
    "P_aero_npeSum;"
    "Events",
    nDeltaBins,
    deltaMin,
    deltaMax,
    nAeroSlices,
    aeroEdges.data()
  );

  auto *hAppliedWeightMap = new TH2D(
    "h_applied_weight_map",
    "Mean applied proton probability;"
    "SHMS #delta [%];"
    "P_aero_npeSum;"
    "#LTw_{p}^{event}#GT",
    nDeltaBins,
    deltaMin,
    deltaMax,
    nAeroSlices,
    aeroEdges.data()
  );

  const std::vector<TH1 *> eventHistograms = {
    hRawMM,
    hProtonMM,
    hCleanedMM,
    hAppliedWeightSumDelta,
    hAppliedWeightCountDelta,
    hAppliedWeightDelta,
    hAppliedWeightSumMap,
    hAppliedWeightCountMap,
    hAppliedWeightMap
  };

  for (
    TH1 *histogram :
    eventHistograms
  ) {
    histogram->SetDirectory(nullptr);
    histogram->Sumw2();
  }

  // ------------------------------------------------------------------
  // Safe event extraction.
  //
  // This replaces the TTreeFormula event loop that caused the crash.
  //
  // V1 = MM
  // V2 = ssdelta
  // V3 = P_aero_npeSum
  // V5 = RF-first branch scan with data-driven RF timing ranges, CTime_ROC1 fallback
  // V8 = single-period RF window selection and corrected RF component ordering
  // V9 = probe RF and CT on every run and select the better validated timing variable
  // ------------------------------------------------------------------

  const Long64_t treeEntries =
    tree->GetEntries();

  tree->SetEstimate(
    treeEntries + 1
  );

  const TString eventExpression =
    TString::Format(
      "%s:%s:%s:%s",
      mmBranch.c_str(),
      deltaBranch.c_str(),
      aeroBranch.c_str(),
      timeBranch.c_str()
    );

  const Long64_t selectedRows =
    tree->Draw(
      eventExpression.Data(),
      analysisCut.c_str(),
      "goff"
    );

  if (selectedRows < 0) {
    std::cerr
      << "TTree::Draw failed while extracting event variables."
      << std::endl;

    inputFile->Close();
    return;
  }

  if (selectedRows == 0) {
    std::cerr
      << "No events passed the event-level selection."
      << std::endl;

    inputFile->Close();
    return;
  }

  const Double_t *mmValues =
    tree->GetV1();

  const Double_t *deltaValues =
    tree->GetV2();

  const Double_t *aeroValues =
    tree->GetV3();

  const Double_t *timingValues =
    tree->GetV4();

  if (
    !mmValues ||
    !deltaValues ||
    !aeroValues ||
    !timingValues
  ) {
    std::cerr
      << "Unable to access TTree::Draw result arrays."
      << std::endl;

    inputFile->Close();
    return;
  }

  Long64_t weightedEvents = 0;
  Long64_t unsupportedEvents = 0;
  Long64_t invalidSliceEvents = 0;

  for (
    Long64_t row = 0;
    row < selectedRows;
    ++row
  ) {
    const double mm =
      mmValues[row];

    const double delta =
      deltaValues[row];

    const double aero =
      aeroValues[row];

    const double timing =
      timingValues[row];

    if (
      !std::isfinite(mm) ||
      !std::isfinite(delta) ||
      !std::isfinite(aero) ||
      !std::isfinite(timing)
    ) {
      continue;
    }

    hRawMM->Fill(
      mm,
      1.0
    );

    int deltaBin =
      static_cast<int>(
        std::floor(
          (delta - deltaMin) /
          deltaBinWidth
        )
      );

    if (delta == deltaMax) {
      deltaBin =
        nDeltaBins - 1;
    }

    if (
      deltaBin < 0 ||
      deltaBin >= nDeltaBins
    ) {
      continue;
    }

    const int aeroSlice =
      findEdgeBin(
        aero,
        aeroEdges
      );

    double protonWeight = 0.0;

    if (
      supportByDelta.at(deltaBin) ==
        SupportClass::Unsupported
    ) {
      ++unsupportedEvents;
    } else if (
      aeroSlice < 0 ||
      aeroSlice >= nAeroSlices ||
      !globalShapes.at(aeroSlice).valid ||
      !deltaSliceFits
         .at(deltaBin)
         .at(aeroSlice)
         .valid
    ) {
      ++invalidSliceEvents;
    } else {
      protonWeight =
        evaluateEventProtonProbability(
          timing,
          globalShapes.at(aeroSlice),
          deltaSliceFits
            .at(deltaBin)
            .at(aeroSlice),
          eventProbabilityDenominatorFloor
        );

      if (protonWeight > 0.0) {
        ++weightedEvents;
      }
    }

    protonWeight =
      std::clamp(
        protonWeight,
        0.0,
        1.0
      );

    hProtonMM->Fill(
      mm,
      protonWeight
    );

    hCleanedMM->Fill(
      mm,
      1.0 - protonWeight
    );

    hAppliedWeightSumDelta->Fill(
      delta,
      protonWeight
    );

    hAppliedWeightCountDelta->Fill(
      delta,
      1.0
    );

    hAppliedWeightSumMap->Fill(
      delta,
      aero,
      protonWeight
    );

    hAppliedWeightCountMap->Fill(
      delta,
      aero,
      1.0
    );
  }

  hAppliedWeightDelta->Divide(
    hAppliedWeightSumDelta,
    hAppliedWeightCountDelta,
    1.0,
    1.0
  );

  hAppliedWeightMap->Divide(
    hAppliedWeightSumMap,
    hAppliedWeightCountMap,
    1.0,
    1.0
  );

  auto *hProtonFractionMM =
    static_cast<TH1D *>(
      hProtonMM->Clone(
        "h_mm_proton_fraction"
      )
    );

  hProtonFractionMM->SetDirectory(nullptr);

  hProtonFractionMM->SetTitle(
    "Estimated proton fraction versus MM;"
    "MM [GeV];"
    "Estimated proton / raw"
  );

  hProtonFractionMM->Divide(hRawMM);

  // ------------------------------------------------------------------
  // Region integrals.
  // ------------------------------------------------------------------

  const int lowMMFirstBin =
    hRawMM
      ->GetXaxis()
      ->FindFixBin(lowMMMin);

  const int lowMMLastBin =
    hRawMM
      ->GetXaxis()
      ->FindFixBin(
        std::nextafter(
          lowMMMax,
          lowMMMin
        )
      );

  const int lambdaFirstBin =
    hRawMM
      ->GetXaxis()
      ->FindFixBin(lambdaMMMin);

  const int lambdaLastBin =
    hRawMM
      ->GetXaxis()
      ->FindFixBin(
        std::nextafter(
          lambdaMMMax,
          lambdaMMMin
        )
      );

  const double rawLowMMIntegral =
    hRawMM->Integral(
      lowMMFirstBin,
      lowMMLastBin
    );

  const double cleanedLowMMIntegral =
    hCleanedMM->Integral(
      lowMMFirstBin,
      lowMMLastBin
    );

  const double rawLambdaIntegral =
    hRawMM->Integral(
      lambdaFirstBin,
      lambdaLastBin
    );

  const double cleanedLambdaIntegral =
    hCleanedMM->Integral(
      lambdaFirstBin,
      lambdaLastBin
    );

  const double lowMMRemovedFraction =
    rawLowMMIntegral > 0.0
      ? 1.0 -
        cleanedLowMMIntegral /
        rawLowMMIntegral
      : 0.0;

  const double lambdaRemovedFraction =
    rawLambdaIntegral > 0.0
      ? 1.0 -
        cleanedLambdaIntegral /
        rawLambdaIntegral
      : 0.0;

  // ------------------------------------------------------------------
  // Global canvas.
  // ------------------------------------------------------------------

  std::vector<TF1 *> componentFunctions;

  auto *globalCanvas = new TCanvas(
    "canvas_global_timing_shapes",
    "Global timing shapes",
    2100,
    1300
  );

  globalCanvas->Divide(3, 2);

  globalCanvas->cd(1);

  gPad->SetRightMargin(0.16);
  gPad->SetLogz();

  hGlobalPID->Draw("COLZ");

  for (
    int aeroSlice = 0;
    aeroSlice < nAeroSlices;
    ++aeroSlice
  ) {
    globalCanvas->cd(
      aeroSlice + 2
    );

    TH1D *projection =
      globalTimingProjections.at(
        aeroSlice
      );

    if (!projection) {
      continue;
    }

    projection->SetLineColor(kBlack);
    projection->SetLineWidth(2);
    projection->Draw("E");

    const TimingShape &shape =
      globalShapes.at(aeroSlice);

    if (shape.fitFunction) {
      shape.fitFunction->SetLineColor(
        shape.valid
          ? kRed + 1
          : kOrange + 7
      );

      shape.fitFunction->SetLineWidth(2);
      shape.fitFunction->Draw("SAME");

      auto *kaonComponent =
        makeGaussianComponent(
          TString::Format(
            "f_global_kaon_component_%d",
            aeroSlice
          ).Data(),
          shape.kaonAmplitude,
          shape.kaonMean,
          shape.kaonSigma,
          timingFitMin,
          timingFitMax,
          kMagenta + 1
        );

      auto *protonComponent =
        makeGaussianComponent(
          TString::Format(
            "f_global_proton_component_%d",
            aeroSlice
          ).Data(),
          shape.protonAmplitude,
          shape.protonMean,
          shape.protonSigma,
          timingFitMin,
          timingFitMax,
          kBlue + 1
        );

      auto *otherComponent =
        makeConstantComponent(
          TString::Format(
            "f_global_other_component_%d",
            aeroSlice
          ).Data(),
          shape.otherAmplitude,
          timingFitMin,
          timingFitMax,
          kGray + 2
        );

      kaonComponent->Draw("SAME");
      protonComponent->Draw("SAME");
      otherComponent->Draw("SAME");

      componentFunctions.push_back(
        kaonComponent
      );

      componentFunctions.push_back(
        protonComponent
      );

      componentFunctions.push_back(
        otherComponent
      );
    }

    auto *label = new TPaveText(
      0.48,
      0.52,
      0.88,
      0.88,
      "NDC"
    );

    label->SetFillStyle(0);
    label->SetBorderSize(0);
    label->SetTextAlign(12);

    label->AddText(
      TString::Format(
        "valid: %s",
        shape.valid ? "yes" : "no"
      )
    );

    label->AddText(
      TString::Format(
        "status: %d",
        shape.fitStatus
      )
    );

    label->AddText(
      TString::Format(
        "bound hit: %s",
        shape.boundHit ? "yes" : "no"
      )
    );

    label->AddText(
      TString::Format(
        "K #mu/#sigma: %.3f / %.3f",
        shape.kaonMean,
        shape.kaonSigma
      )
    );

    label->AddText(
      TString::Format(
        "p #mu/#sigma: %.3f / %.3f",
        shape.protonMean,
        shape.protonSigma
      )
    );

    label->AddText(
      TString::Format(
        "separation: %.2f",
        shape.separation
      )
    );

    label->AddText(
      TString::Format(
        "K/p significance: %.1f / %.1f",
        shape.kaonSignificance,
        shape.protonSignificance
      )
    );

    label->AddText(
      TString::Format(
        "Poisson D/ndf: %.2f (%d)",
        shape.poissonDevianceNdf,
        shape.goodnessNdf
      )
    );

    label->AddText(
      TString::Format(
        "D/N: %.3f",
        shape.poissonDeviancePerEntry
      )
    );

    label->Draw();
  }

  // ------------------------------------------------------------------
  // Summary canvas.
  // ------------------------------------------------------------------

  auto *summaryCanvas = new TCanvas(
    "canvas_proton_cleaning_summary",
    "Event-level proton-cleaning summary",
    2100,
    1300
  );

  summaryCanvas->Divide(3, 2);

  summaryCanvas->cd(1);

  hFitProtonWeightDelta->SetMinimum(0.0);
  hFitProtonWeightDelta->SetMaximum(1.05);

  hFitProtonWeightDelta->SetLineColor(
    kGray + 2
  );

  hFitProtonWeightDelta->SetMarkerColor(
    kGray + 2
  );

  hFitProtonWeightDelta->SetMarkerStyle(20);

  hAppliedWeightDelta->SetLineColor(
    kRed + 1
  );

  hAppliedWeightDelta->SetMarkerColor(
    kRed + 1
  );

  hAppliedWeightDelta->SetMarkerStyle(21);
  hAppliedWeightDelta->SetLineWidth(3);

  hFitProtonWeightDelta->Draw("E1");
  hAppliedWeightDelta->Draw("E1 SAME");

  auto *weightLegend = new TLegend(
    0.47,
    0.72,
    0.88,
    0.88
  );

  weightLegend->AddEntry(
    hFitProtonWeightDelta,
    "integrated fitted fraction",
    "lep"
  );

  weightLegend->AddEntry(
    hAppliedWeightDelta,
    "mean applied event weight",
    "lep"
  );

  weightLegend->Draw();

  summaryCanvas->cd(2);

  hProtonYield->SetLineColor(kRed + 1);
  hProtonYield->SetLineWidth(3);

  hKaonYield->SetLineColor(kMagenta + 1);
  hKaonYield->SetLineWidth(3);

  hOtherYield->SetLineColor(kBlue + 1);
  hOtherYield->SetLineWidth(3);

  const double maximumYield =
    std::max(
      {
        hProtonYield->GetMaximum(),
        hKaonYield->GetMaximum(),
        hOtherYield->GetMaximum()
      }
    );

  hProtonYield->SetMaximum(
    1.25 * maximumYield
  );

  hProtonYield->Draw("HIST");
  hKaonYield->Draw("HIST SAME");
  hOtherYield->Draw("HIST SAME");

  auto *yieldLegend = new TLegend(
    0.57,
    0.70,
    0.88,
    0.88
  );

  yieldLegend->AddEntry(
    hProtonYield,
    "proton",
    "l"
  );

  yieldLegend->AddEntry(
    hKaonYield,
    "kaon",
    "l"
  );

  yieldLegend->AddEntry(
    hOtherYield,
    "other",
    "l"
  );

  yieldLegend->Draw();

  summaryCanvas->cd(3);

  gPad->SetRightMargin(0.16);

  hAppliedWeightMap->SetMinimum(0.0);
  hAppliedWeightMap->SetMaximum(1.0);
  hAppliedWeightMap->Draw("COLZ TEXT");

  summaryCanvas->cd(4);

  hFitCoverage->SetLineColor(kBlue + 1);
  hFitCoverage->SetMarkerColor(kBlue + 1);
  hFitCoverage->SetMarkerStyle(20);
  hFitCoverage->SetMinimum(0.0);
  hFitCoverage->SetMaximum(1.05);
  hFitCoverage->Draw("E1");

  auto *supportedCoverageLine =
    new TLine(
      deltaMin,
      minimumSupportedCoverage,
      deltaMax,
      minimumSupportedCoverage
    );

  supportedCoverageLine->SetLineStyle(2);
  supportedCoverageLine->SetLineColor(
    kGreen + 2
  );

  supportedCoverageLine->Draw("SAME");

  auto *marginalCoverageLine =
    new TLine(
      deltaMin,
      minimumMarginalCoverage,
      deltaMax,
      minimumMarginalCoverage
    );

  marginalCoverageLine->SetLineStyle(2);
  marginalCoverageLine->SetLineColor(
    kOrange + 7
  );

  marginalCoverageLine->Draw("SAME");

  summaryCanvas->cd(5);

  hRawMM->SetLineColor(kBlack);
  hRawMM->SetLineWidth(3);

  hProtonMM->SetLineColor(kRed + 1);
  hProtonMM->SetLineWidth(3);

  hCleanedMM->SetLineColor(kGreen + 2);
  hCleanedMM->SetLineWidth(3);

  const double maximumMM =
    std::max(
      {
        hRawMM->GetMaximum(),
        hProtonMM->GetMaximum(),
        hCleanedMM->GetMaximum()
      }
    );

  hRawMM->SetMaximum(
    1.20 * maximumMM
  );

  hRawMM->Draw("HIST");
  hProtonMM->Draw("HIST SAME");
  hCleanedMM->Draw("HIST SAME");

  auto *mmLegend = new TLegend(
    0.47,
    0.68,
    0.88,
    0.88
  );

  mmLegend->AddEntry(
    hRawMM,
    "raw kaon-selected MM",
    "l"
  );

  mmLegend->AddEntry(
    hProtonMM,
    "estimated proton contamination",
    "l"
  );

  mmLegend->AddEntry(
    hCleanedMM,
    "proton-cleaned kaon MM",
    "l"
  );

  mmLegend->Draw();

  auto *integralText = new TPaveText(
    0.48,
    0.47,
    0.88,
    0.65,
    "NDC"
  );

  integralText->SetFillColor(kWhite);
  integralText->SetFillStyle(1001);
  integralText->SetBorderSize(1);
  integralText->SetTextAlign(12);

  integralText->AddText(
    TString::Format(
      "Low-MM removed: %.1f%%",
      100.0 * lowMMRemovedFraction
    )
  );

  integralText->AddText(
    TString::Format(
      "#Lambda region removed: %.1f%%",
      100.0 * lambdaRemovedFraction
    )
  );

  integralText->Draw();

  summaryCanvas->cd(6);

  hProtonFractionMM->SetLineColor(kRed + 1);
  hProtonFractionMM->SetLineWidth(3);
  hProtonFractionMM->SetMinimum(0.0);
  hProtonFractionMM->SetMaximum(1.05);
  hProtonFractionMM->Draw("HIST");

  // ------------------------------------------------------------------
  // Per-delta canvases.
  // ------------------------------------------------------------------

  std::vector<TCanvas *> deltaCanvases(
    nDeltaBins,
    nullptr
  );

  for (
    int deltaBin = 0;
    deltaBin < nDeltaBins;
    ++deltaBin
  ) {
    const double deltaLow =
      deltaMin +
      deltaBin * deltaBinWidth;

    const double deltaHigh =
      deltaLow + deltaBinWidth;

    auto *canvas = new TCanvas(
      TString::Format(
        "canvas_delta_fit_%d",
        deltaBin
      ),
      TString::Format(
        "Delta-bin fit %d",
        deltaBin
      ),
      2100,
      1300
    );

    canvas->Divide(3, 2);

    canvas->cd(1);

    gPad->SetRightMargin(0.16);

    if (deltaPIDHistograms.at(deltaBin)) {
      gPad->SetLogz();

      deltaPIDHistograms
        .at(deltaBin)
        ->Draw("COLZ");
    }

    auto *summaryText = new TPaveText(
      0.43,
      0.53,
      0.88,
      0.88,
      "NDC"
    );

    summaryText->SetFillColor(kWhite);
    summaryText->SetFillStyle(1001);
    summaryText->SetBorderSize(2);
    summaryText->SetTextAlign(12);

    summaryText->SetTextColor(
      supportClassColor(
        supportByDelta.at(deltaBin)
      )
    );

    summaryText->AddText(
      TString::Format(
        "%.1f #leq #delta < %.1f",
        deltaLow,
        deltaHigh
      )
    );

    summaryText->AddText(
      TString::Format(
        "support: %s",
        supportClassLabel(
          supportByDelta.at(deltaBin)
        )
      )
    );

    summaryText->AddText(
      TString::Format(
        "N_{p} = %.1f",
        protonYieldByDelta.at(deltaBin)
      )
    );

    summaryText->AddText(
      TString::Format(
        "N_{K} = %.1f",
        kaonYieldByDelta.at(deltaBin)
      )
    );

    summaryText->AddText(
      TString::Format(
        "N_{other} = %.1f",
        otherYieldByDelta.at(deltaBin)
      )
    );

    summaryText->AddText(
      TString::Format(
        "fit coverage = %.3f",
        validCoverageByDelta.at(deltaBin)
      )
    );

    summaryText->AddText(
      TString::Format(
        "valid slices = %d/%d",
        validSlicesByDelta.at(deltaBin),
        nAeroSlices
      )
    );

    summaryText->AddText(
      TString::Format(
        "fit-integrated w_{p} = %.3f",
        hFitProtonWeightDelta
          ->GetBinContent(deltaBin + 1)
      )
    );

    summaryText->Draw();

    for (
      int aeroSlice = 0;
      aeroSlice < nAeroSlices;
      ++aeroSlice
    ) {
      canvas->cd(
        aeroSlice + 2
      );

      TH1D *projection =
        deltaTimingProjections
          .at(deltaBin)
          .at(aeroSlice);

      if (!projection) {
        continue;
      }

      projection->SetLineColor(kBlack);
      projection->SetLineWidth(2);
      projection->Draw("E");

      const SliceFitResult &sliceFit =
        deltaSliceFits
          .at(deltaBin)
          .at(aeroSlice);

      const TimingShape &shape =
        globalShapes.at(aeroSlice);

      if (sliceFit.fitFunction) {
        sliceFit.fitFunction->SetLineColor(
          sliceFit.valid
            ? kRed + 1
            : kOrange + 7
        );

        sliceFit.fitFunction->SetLineWidth(2);
        sliceFit.fitFunction->Draw("SAME");

        auto *kaonComponent =
          makeGaussianComponent(
            TString::Format(
              "f_delta_%d_aero_%d_kaon",
              deltaBin,
              aeroSlice
            ).Data(),
            sliceFit.kaonAmplitude,
            shape.kaonMean,
            shape.kaonSigma,
            timingFitMin,
            timingFitMax,
            kMagenta + 1
          );

        auto *protonComponent =
          makeGaussianComponent(
            TString::Format(
              "f_delta_%d_aero_%d_proton",
              deltaBin,
              aeroSlice
            ).Data(),
            sliceFit.protonAmplitude,
            shape.protonMean,
            shape.protonSigma,
            timingFitMin,
            timingFitMax,
            kBlue + 1
          );

        auto *otherComponent =
          makeConstantComponent(
            TString::Format(
              "f_delta_%d_aero_%d_other",
              deltaBin,
              aeroSlice
            ).Data(),
            sliceFit.otherAmplitude,
            timingFitMin,
            timingFitMax,
            kGray + 2
          );

        kaonComponent->Draw("SAME");
        protonComponent->Draw("SAME");
        otherComponent->Draw("SAME");

        componentFunctions.push_back(
          kaonComponent
        );

        componentFunctions.push_back(
          protonComponent
        );

        componentFunctions.push_back(
          otherComponent
        );
      }

      auto *fitText = new TPaveText(
        0.49,
        0.52,
        0.88,
        0.88,
        "NDC"
      );

      fitText->SetFillStyle(0);
      fitText->SetBorderSize(0);
      fitText->SetTextAlign(12);

      fitText->AddText(
        TString::Format(
          "global shape valid: %s",
          shape.valid ? "yes" : "no"
        )
      );

      fitText->AddText(
        TString::Format(
          "slice fit valid: %s",
          sliceFit.valid ? "yes" : "no"
        )
      );

      fitText->AddText(
        TString::Format(
          "status: %d",
          sliceFit.fitStatus
        )
      );

      fitText->AddText(
        TString::Format(
          "N_{p}: %.1f",
          sliceFit.protonYield
        )
      );

      fitText->AddText(
        TString::Format(
          "N_{K}: %.1f",
          sliceFit.kaonYield
        )
      );

      fitText->AddText(
        TString::Format(
          "N_{other}: %.1f",
          sliceFit.otherYield
        )
      );

      fitText->AddText(
        TString::Format(
          "model/data: %.2f",
          sliceFit.modelDataRatio
        )
      );

      fitText->AddText(
        TString::Format(
          "Poisson D/ndf: %.2f",
          sliceFit.poissonDevianceNdf
        )
      );

      fitText->AddText(
        TString::Format(
          "D/N: %.3f",
          sliceFit.poissonDeviancePerEntry
        )
      );

      fitText->Draw();
    }

    deltaCanvases.at(deltaBin) =
      canvas;
  }

  // ------------------------------------------------------------------
  // One multipage PDF.
  // ------------------------------------------------------------------

  globalCanvas->Modified();
  globalCanvas->Update();

  summaryCanvas->Modified();
  summaryCanvas->Update();

  globalCanvas->Print(
    (outputPDF + "[").c_str()
  );

  globalCanvas->Print(
    outputPDF.c_str()
  );

  summaryCanvas->Print(
    outputPDF.c_str()
  );

  for (
    TCanvas *canvas :
    deltaCanvases
  ) {
    if (!canvas) {
      continue;
    }

    canvas->Modified();
    canvas->Update();

    canvas->Print(
      outputPDF.c_str()
    );
  }

  summaryCanvas->Print(
    (outputPDF + "]").c_str()
  );

  // ------------------------------------------------------------------
  // One ROOT file.
  // ------------------------------------------------------------------

  TFile outputFile(
    outputROOT.c_str(),
    "RECREATE"
  );

  if (outputFile.IsZombie()) {
    std::cerr
      << "Unable to create output file: "
      << outputROOT
      << std::endl;

    inputFile->Close();
    return;
  }

  TDirectory *globalDirectory =
    outputFile.mkdir("global");

  globalDirectory->cd();

  hGlobalPID->Write();

  for (
    int aeroSlice = 0;
    aeroSlice < nAeroSlices;
    ++aeroSlice
  ) {
    if (globalTimingProjections.at(aeroSlice)) {
      globalTimingProjections
        .at(aeroSlice)
        ->Write();
    }

    if (globalShapes.at(aeroSlice).fitFunction) {
      globalShapes
        .at(aeroSlice)
        .fitFunction
        ->Write();
    }
  }

  globalCanvas->Write();

  outputFile.cd();

  TDirectory *summaryDirectory =
    outputFile.mkdir("summary");

  summaryDirectory->cd();

  hProtonYield->Write();
  hKaonYield->Write();
  hOtherYield->Write();

  hFitProtonWeightDelta->Write();
  hAppliedWeightDelta->Write();
  hAppliedWeightMap->Write();

  hFitChi2->Write();
  hFitCoverage->Write();
  hValidSlices->Write();
  hSupportClass->Write();

  hRawMM->Write();
  hProtonMM->Write();
  hCleanedMM->Write();
  hProtonFractionMM->Write();

  summaryCanvas->Write();

  outputFile.cd();

  TDirectory *deltaDirectory =
    outputFile.mkdir("delta_bins");

  for (
    int deltaBin = 0;
    deltaBin < nDeltaBins;
    ++deltaBin
  ) {
    deltaDirectory->cd();

    TDirectory *singleDeltaDirectory =
      deltaDirectory->mkdir(
        TString::Format(
          "delta_bin_%02d",
          deltaBin
        ).Data()
      );

    singleDeltaDirectory->cd();

    TParameter<int> supportParameter(
      "support_class",
      static_cast<int>(
        supportByDelta.at(deltaBin)
      )
    );

    supportParameter.Write();

    TNamed supportLabel(
      "support_label",
      supportClassLabel(
        supportByDelta.at(deltaBin)
      )
    );

    supportLabel.Write();

    if (deltaPIDHistograms.at(deltaBin)) {
      deltaPIDHistograms
        .at(deltaBin)
        ->Write();
    }

    for (
      int aeroSlice = 0;
      aeroSlice < nAeroSlices;
      ++aeroSlice
    ) {
      TH1D *projection =
        deltaTimingProjections
          .at(deltaBin)
          .at(aeroSlice);

      if (projection) {
        projection->Write();
      }

      TF1 *fitFunction =
        deltaSliceFits
          .at(deltaBin)
          .at(aeroSlice)
          .fitFunction;

      if (fitFunction) {
        fitFunction->Write();
      }
    }

    if (deltaCanvases.at(deltaBin)) {
      deltaCanvases
        .at(deltaBin)
        ->Write();
    }
  }

  outputFile.cd();

  TDirectory *metadataDirectory =
    outputFile.mkdir("metadata");

  metadataDirectory->cd();

  TNamed macroVersionTag(
    "macro_version",
    macroVersion
  );
  macroVersionTag.Write();

  TNamed timingBranchTag(
    "timing_branch",
    timeBranch.c_str()
  );
  timingBranchTag.Write();

  TNamed rfTimingBranchTag(
    "rf_timing_branch",
    rfTimeBranch.c_str()
  );
  rfTimingBranchTag.Write();

  TParameter<int>(
    "rf_timing_attempted",
    rfTimingAttempted ? 1 : 0
  ).Write();

  TParameter<int>(
    "rf_timing_selected",
    rfTimingSelected ? 1 : 0
  ).Write();

  TParameter<int>(
    "proton_peak_is_lower",
    protonPeakIsLower ? 1 : 0
  ).Write();

  TParameter<int>(
    "ct_fallback_used",
    ctFallbackUsed ? 1 : 0
  ).Write();

  TParameter<int>(
    "ct_timing_evaluated",
    ctTimingEvaluated ? 1 : 0
  ).Write();

  TParameter<int>(
    "rf_probe_valid_global_shapes",
    rfProbeValidShapes
  ).Write();

  TParameter<int>(
    "ct_probe_valid_global_shapes",
    ctProbeValidShapes
  ).Write();

  TNamed timingSelectionReasonTag(
    "timing_selection_reason",
    timingSelectionReason.c_str()
  );
  timingSelectionReasonTag.Write();

  TNamed fitStatisticTag(
    "fit_statistic",
    "Baker-Cousins Poisson deviance; bin-integrated likelihood"
  );
  fitStatisticTag.Write();

  TNamed fitValidationStatisticTag(
    "fit_validation_statistic",
    useDeviancePerEntryValidation
      ? "poisson_deviance_per_entry"
      : "poisson_deviance_per_ndf"
  );
  fitValidationStatisticTag.Write();

  TParameter<double>(
    "maximum_global_poisson_deviance_per_ndf",
    maximumGlobalPoissonDevianceNdf
  ).Write();

  TParameter<double>(
    "maximum_global_poisson_deviance_per_entry",
    maximumGlobalPoissonDeviancePerEntry
  ).Write();

  TParameter<double>(
    "maximum_slice_poisson_deviance_per_ndf",
    maximumSlicePoissonDevianceNdf
  ).Write();

  TParameter<double>(
    "maximum_slice_poisson_deviance_per_entry",
    maximumSlicePoissonDeviancePerEntry
  ).Write();

  TParameter<int>(
    "rf_correction_used",
    0
  ).Write();

  TParameter<double>(
    "beam_bunch_spacing_ns",
    beamBunchSpacingNs
  ).Write();

  TParameter<double>(
    "timing_fit_min_ns",
    timingFitMin
  ).Write();

  TParameter<double>(
    "timing_fit_max_ns",
    timingFitMax
  ).Write();

  TParameter<Long64_t>(
    "selected_events",
    selectedRows
  ).Write();

  TParameter<Long64_t>(
    "weighted_events",
    weightedEvents
  ).Write();

  TParameter<Long64_t>(
    "unsupported_delta_events",
    unsupportedEvents
  ).Write();

  TParameter<Long64_t>(
    "invalid_pid_slice_events",
    invalidSliceEvents
  ).Write();

  TParameter<double>(
    "low_mm_removed_fraction",
    lowMMRemovedFraction
  ).Write();

  TParameter<double>(
    "lambda_removed_fraction",
    lambdaRemovedFraction
  ).Write();

  outputFile.Write();
  outputFile.Close();

  inputFile->Close();

  std::cout
    << "\n===== Event-level proton-cleaning POC ====="
    << "\nInput file: "
    << filename
    << "\nTree: "
    << treeName
    << "\nTiming selection: "
    << timeBranch
    << " ("
    << timingSelectionReason
    << ")"
    << "\nGlobal PID entries: "
    << hGlobalPID->Integral()
    << "\nIdentifiable global aerogel slices: "
    << validGlobalShapes
    << " / "
    << nAeroSlices
    << "\nSelected events: "
    << selectedRows
    << "\nEvents receiving nonzero proton weight: "
    << weightedEvents
    << "\nEvents in unsupported delta bins: "
    << unsupportedEvents
    << "\nEvents in invalid PID slices: "
    << invalidSliceEvents
    << "\nRaw MM integral: "
    << hRawMM->Integral()
    << "\nEstimated proton MM integral: "
    << hProtonMM->Integral()
    << "\nProton-cleaned MM integral: "
    << hCleanedMM->Integral()
    << "\nLow-MM removed fraction: "
    << 100.0 * lowMMRemovedFraction
    << "%"
    << "\nLambda-region removed fraction: "
    << 100.0 * lambdaRemovedFraction
    << "%"
    << "\n\nCreated exactly two files:"
    << "\n  PDF:  "
    << outputPDF
    << "\n  ROOT: "
    << outputROOT
    << "\n==========================================="
    << std::endl;
}
