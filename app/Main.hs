{-# LANGUAGE NoFieldSelectors #-}
{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE StrictData #-}
{-# LANGUAGE ViewPatterns #-}
{-# LANGUAGE LambdaCase #-}
module Main where

import Control.Monad
import Data.List (intercalate)
import Data.Number.Erf
import Data.Random
import Data.Random.Distribution.Binomial (binomial)
import Statistics.Distribution (quantile, complQuantile)
import Statistics.Distribution.Beta (betaDistr)
import System.Random.Stateful
import System.Environment

type Formula = Double -> Int -> Int -> (Double, Double)

dumbInterval :: Formula
dumbInterval _ _ _ = (0, 1)

waldInterval :: Formula
waldInterval z (fromIntegral -> p) (fromIntegral -> n) = (phat - term, phat + term)
  where
    phat = p / n
    term = z * sqrt (phat * (1 - phat) / n)

wilsonScoreInterval :: Formula
wilsonScoreInterval z sI nI
  | sI == 0 = (0, upper)
  | sI == nI = (lower, 1)
  | otherwise = (lower, upper)
  where
    s = fromIntegral sI
    n = fromIntegral nI
    z2 = z * z
    phat = s / n
    term = z / (2 * n) * sqrt (4 * n * phat * (1 - phat) + z2)
    factor = 1 / (1 + z2 / n)
    offset = z2 / (2 * n)
    lower = factor * (phat + offset - term)
    upper = factor * (phat + offset + term)

wilsonScoreIntervalCC :: Formula
wilsonScoreIntervalCC z s num
  | s == 0 = (0, upper)
  | s == num = (lower, 1)
  | otherwise = (lower, upper)
  where
    n = fromIntegral num
    z2 = z * z
    phat = fromIntegral s / n
    a = 2 * n * phat + z2
    t = z2 - 1 / n + 4 * n * phat * (1 - phat)
    v = 4 * phat - 2
    bL = z * sqrt (t + v) + 1
    bU = z * sqrt (t - v) + 1
    c = 2 * (n + z2)
    lower = max 0 $ (a - bL) / c
    upper = min 1 $ (a + bU) / c

agrestiCoullInterval :: Formula
agrestiCoullInterval z (fromIntegral -> s) (fromIntegral -> n) = (p' - d, p' + d)
  where
    n' = n + z * z
    p' = (s + z * z / 2) / n'
    d = z * sqrt (p' * (1 - p') / n')

clopperPearsonInterval :: Formula
clopperPearsonInterval z sI nI
  | sI == 0 = (0, 1 - alpha ** (1 / n))
  | sI == nI = (alpha ** (1 / n), 1)
  | otherwise = (quantile distrL alpha, complQuantile distrU alpha)
  where
    s = fromIntegral sI
    n = fromIntegral nI
    alpha = erfc (z / sqrt 2) / 2
    distrL = betaDistr s (n - s + 1)
    distrU = betaDistr (s + 1) (n - s)

data Session = Session
  { numTrials :: Int
  , numSuccess :: Int
  , lower :: Double
  , upper :: Double
  }

sessionR :: Formula -> Double -> Double -> RVar Session
sessionR func zScore prob = do
  numTrials <- Data.Random.uniform 1 1000
  numSuccess <- binomial (numTrials :: Int) prob
  let (lower, upper) = func zScore numSuccess numTrials
  pure Session {
    numTrials,
    numSuccess,
    lower,
    upper
  }

data Evaluation = Evaluation
  { rate :: Double
  , score :: Double
  }

evaluateR :: Formula -> Double -> Double -> RVar Evaluation
evaluateR func zScore prob = do
  let numSessions = 10000
  intervals <- replicateM numSessions $ sessionR func zScore prob
  let count = sum [1 :: Int | Session{lower, upper} <- intervals, lower <= prob, prob <= upper]
  let rate = fromIntegral count / fromIntegral numSessions
  let sumWidth = sum [fromIntegral numTrials * (upper - lower) ^ (2 :: Int) | Session{..} <- intervals]
  pure Evaluation { rate, score = sumWidth / fromIntegral numSessions }

significanceToZ :: Double -> Double
significanceToZ sig = sqrt 2 * inverfc sig

algorithms :: [(String, Formula)]
algorithms = 
  [ 
  -- ("dumb", dumbInterval)
    ("clopper-pearson", clopperPearsonInterval)
  , ("wald", waldInterval)
  , ("wilson-score", wilsonScoreInterval)
  , ("wilson-score-cc", wilsonScoreIntervalCC)
  , ("agresti-coull", agrestiCoullInterval)
  ]

runR :: RVar a -> IO a
runR rv = initStdGen >>= newIOGenM >>= runRVar rv

main :: IO ()
main = getArgs >>= \case
  ["samples"] -> do
    putStrLn "algorithm,prob,numTrials,numSuccess,width"
    replicateM_ 1000 $ do
      forM_ algorithms $ \(name, func) -> do
        (prob, Session{..}) <- runR $ do
          prob <- Data.Random.uniform 0 1
          session <- sessionR func 1.96 prob
          pure (prob, session)
        putStrLn $ intercalate "," [name, show prob, show numTrials, show numSuccess, show $ upper - lower]
  ["evaluate"] -> do
    putStrLn "algorithm,prob,rate,avgWidth"
    forM_ algorithms $ \(name, func) -> do
      forM_ [0,0.01..1] $ \prob -> do
        gen <- newIOGenM =<< initStdGen
        let significance = 0.05
        let zScore = significanceToZ significance
        Evaluation { rate,score } <- flip runRVar gen $ evaluateR func zScore prob
        putStrLn $ intercalate "," [name, show prob, show rate, show score]
  _ -> putStrLn "Usage: binomial-confidence-interval samples|evaluate"
