module Mandelbrot.Graphics.Trustworthy.Quadratic
  ( quadratic
  ) where

import qualified Mandelbrot.Graphics.BoundingBox.D1 as D1
import Mandelbrot.Graphics.BoundingBox.D2

quadratic :: (Real r) => Box r -> Box r -> Box r
{-# SPECIALIZE quadratic :: Box Double -> Box Double -> Box Double #-}
{-# SPECIALIZE quadratic :: Box Rational -> Box Rational -> Box Rational #-}
quadratic (Box (D1.Box alo ahi) (D1.Box blo bhi)) (Box (D1.Box xlo xhi) (D1.Box ylo yhi)) =
    Box (D1.Box (minimum us) (maximum us)) (D1.Box (minimum vs) (maximum vs))
  where
    (us, vs) = unzip
      [ (x^2-y^2+a, 2*x*y+b)
      | x <- [xlo, xhi]
      , y <- [ylo, yhi]
      , a <- [alo, ahi]
      , b <- [blo, bhi]
      ]
