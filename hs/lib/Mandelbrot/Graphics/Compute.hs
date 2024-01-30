module Mandelbrot.Graphics.Compute
  ( Compute(..)
  , Partial(..)
  , compute
  , initialize
  , step
  , fixup
  , isInterior
  , interiorDE
  ) where

import Data.Complex (Complex())
import Data.Maybe (isNothing, isJust)

import Mandelbrot.Numerics (Square, Approx, attractor_, lastMay, finite, (!), magnitude', magnitudeSquared)

data Partial r = Partial !Int !(Complex r)
  deriving Show

data Compute r
  = Interior{ computeP :: !Int, computeDZ :: !(Complex r), computeDE :: !r }
  | Exterior{ computeP, computeN :: !Int, computeZ :: !(Complex r), computeDE :: !r }
  | Unknown{ computePartials :: Maybe [Partial r], computeP, computeN :: !Int, computeC, computeDC0, computeZ, computeDC :: !(Complex r), computeER2, computeMZ2 :: !r }
  deriving Show

initialize :: RealFloat r => Bool -> r -> (Complex r, Complex r) -> Compute r
initialize bias er (c, dc) = Unknown
  { computePartials = if bias then Nothing else Just []
  , computeP = 0
  , computeN = 0
  , computeC = c
  , computeDC0 = dc
  , computeZ = 0
  , computeDC = 0
  , computeER2 = er * er
  , computeMZ2 = 1 / 0
  }

step :: (RealFloat r, Square r, Approx r) => Compute r -> Compute r
step u@Unknown{ computePartials = mps, computeN = n, computeC = c, computeDC0 = dc00, computeZ = z0, computeDC = dc0, computeER2 = er2, computeMZ2 = mz2 }
  | z2 < mz2 / 4 && isNothing mps && isJust interior = case interior of
      Just (de, dz) -> Interior{ computeP = p, computeDZ = dz, computeDE = de / magnitude' dc00 }
      Nothing -> error "step impossible"
  | z2 < mz2 && z2 < er2 = u{ computePartials = (if z2 < mz2 / 4 then (Partial p z :) else id) `fmap` mps, computeP = p, computeN = p, computeZ = z, computeDC = dc, computeMZ2 = z2 }
  | z2 < er2 = u{ computeN = p, computeZ = z, computeDC = dc }
  | otherwise = Exterior{ computeP = computeP u, computeN = p, computeZ = z, computeDE = exterior }
  where
    dc = 2 * z0 * dc0 + dc00
    z = z0 * z0 + c
    z2 = magnitudeSquared z
    p = n + 1
    interior = interiorDE z c p 64
    exterior = sqrt z2 * log z2 / magnitude' dc
step u = u

fixup :: (RealFloat r, Square r, Approx r) => Compute r -> Compute r
fixup u@Unknown{ computePartials = Just ps, computeC = c } = go (reverse ps)
  where
    go [] = u{ computePartials = Just [] }
    go (Partial p z : pzs) = case interiorDE z c p 64 of
      Just (de, dz) -> Interior{ computeP = p, computeDZ = dz, computeDE = de }
      Nothing -> go pzs
fixup u = u

compute :: (RealFloat r, Square r, Approx r) => Bool -> r -> (Complex r, Complex r) -> Int -> Compute r
compute bias er cdc = go (initialize bias er cdc)
  where
    go r@Exterior{} _ = r
    go r@Interior{} _ = r
    go r 0 = fixup r
    go r n = go (step r) (n - 1)

isInterior :: Compute r -> Bool
isInterior Interior{} = True
isInterior _ = False

interiorDE :: (RealFloat r, Square r, Approx r) => Complex r -> Complex r -> Int -> Int -> Maybe (r, Complex r)
interiorDE z0 c p steps = case lastMay . take steps $ attractor_ p c z0 of
  Just z00 -> case go1 p z00 1 of
    dz0 | finite dz0 && magnitudeSquared dz0 <= 1 -> Just (go2 p z00 1 0 0 0)
    _ -> Nothing
  _ -> Nothing
  where
    go1 q z dz
      | q > 0 = go1 ! (q - 1) ! (z * z + c) ! (2 * z * dz)
      | otherwise = dz
    go2 q z dz dzdz dc dcdz
      | q > 0 = go2 ! (q - 1) ! (z * z + c) ! (2 * z * dz) ! (2 * (dz * dz + z * dzdz)) ! (2 * z * dc + 1) ! (2 * (z * dcdz + dz * dc))
      | otherwise = ((1 - magnitudeSquared dz) / magnitude' (dcdz + dzdz * dc / (1 - dz)), dz)
