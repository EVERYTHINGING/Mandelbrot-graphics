{-# LANGUAGE DeriveFunctor #-}
{-# LANGUAGE DeriveFoldable #-}
{-# LANGUAGE DeriveTraversable #-}
module Mandelbrot.Graphics.BoundingBox.D2 where

import qualified Mandelbrot.Graphics.BoundingBox.D1 as D1

data Box r = Box
  { xRange :: !(D1.Box r)
  , yRange :: !(D1.Box r)
  }
  deriving (Eq, Read, Show, Functor, Foldable, Traversable)

split :: RealFrac r => Box r -> (Box r, Box r, Box r, Box r)
{-# SPECIALIZE split :: Box Double -> (Box Double, Box Double, Box Double, Box Double) #-}
{-# SPECIALIZE split :: Box Rational -> (Box Rational, Box Rational, Box Rational, Box Rational) #-}
split (Box x y) = (Box xlo ylo, Box xlo yhi, Box xhi ylo, Box xhi yhi)
  where
    (xlo, xhi) = D1.split x
    (ylo, yhi) = D1.split y

inside  :: Ord r => Box r -> Box r -> Bool
{-# SPECIALIZE inside :: Box Double -> Box Double -> Bool #-}
{-# SPECIALIZE inside :: Box Rational -> Box Rational -> Bool #-}
inside  (Box x y) (Box u v) = x `D1.inside`  u && y `D1.inside`  v

inside' :: Ord r => Box r -> Box r -> Bool
{-# SPECIALIZE inside' :: Box Double -> Box Double -> Bool #-}
{-# SPECIALIZE inside' :: Box Rational -> Box Rational -> Bool #-}
inside' (Box x y) (Box u v) = x `D1.inside'` u && y `D1.inside'` v

disjoint  :: Ord r => Box r -> Box r -> Bool
{-# SPECIALIZE disjoint :: Box Double -> Box Double -> Bool #-}
{-# SPECIALIZE disjoint :: Box Rational -> Box Rational -> Bool #-}
disjoint  (Box x y) (Box u v) = x `D1.disjoint`  u || y `D1.disjoint`  v

disjoint' :: Ord r => Box r -> Box r -> Bool
{-# SPECIALIZE disjoint' :: Box Double -> Box Double -> Bool #-}
{-# SPECIALIZE disjoint' :: Box Rational -> Box Rational -> Bool #-}
disjoint' (Box x y) (Box u v) = x `D1.disjoint'` u || y `D1.disjoint'` v

symmetric :: Num r => r -> Box r
{-# SPECIALIZE symmetric :: Double -> Box Double #-}
{-# SPECIALIZE symmetric :: Rational -> Box Rational #-}
symmetric r = Box b b where b = D1.symmetric r

singleton :: r -> r -> Box r
{-# SPECIALIZE singleton :: Double -> Double -> Box Double #-}
{-# SPECIALIZE singleton :: Rational -> Rational -> Box Rational #-}
singleton x y = Box (D1.singleton x) (D1.singleton y)

union :: Ord r => Box r -> Box r -> Box r
{-# SPECIALIZE union :: Box Double -> Box Double -> Box Double #-}
{-# SPECIALIZE union :: Box Rational -> Box Rational -> Box Rational #-}
union (Box a b) (Box u v) = Box (D1.union a u) (D1.union b v)
