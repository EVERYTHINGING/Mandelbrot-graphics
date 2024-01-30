{-# LANGUAGE DeriveFunctor #-}
{-# LANGUAGE DeriveFoldable #-}
{-# LANGUAGE DeriveTraversable #-}
module Mandelbrot.Graphics.BoundingBox.D1 where

data Box r = Box
  { lower :: !r
  , upper :: !r
  }
  deriving (Eq, Read, Show, Functor, Foldable, Traversable)

box :: Ord r => r -> r -> Box r
{-# SPECIALIZE box :: Double -> Double -> Box Double #-}
{-# SPECIALIZE box :: Rational -> Rational -> Box Rational #-}
box a b
  | a <= b    = Box a b
  | otherwise = Box b a

midPoint :: Fractional r => Box r -> r
{-# SPECIALIZE midPoint :: Box Double -> Double #-}
{-# SPECIALIZE midPoint :: Box Rational -> Rational #-}
midPoint (Box a b) = (a + b) / 2

split :: RealFrac r => Box r -> (Box r, Box r)
{-# SPECIALIZE split :: Box Double -> (Box Double, Box Double) #-}
{-# SPECIALIZE split :: Box Rational -> (Box Rational, Box Rational) #-}
split x@(Box a b) = (Box a m, Box m b) where m = midPoint x

inside  :: Ord r => Box r -> Box r -> Bool
{-# SPECIALIZE inside :: Box Double -> Box Double -> Bool #-}
{-# SPECIALIZE inside :: Box Rational -> Box Rational -> Bool #-}
inside  (Box xlo xhi) (Box ulo uhi) = ulo <= xlo && xhi <= uhi

inside' :: Ord r => Box r -> Box r -> Bool
{-# SPECIALIZE inside' :: Box Double -> Box Double -> Bool #-}
{-# SPECIALIZE inside' :: Box Rational -> Box Rational -> Bool #-}
inside' (Box xlo xhi) (Box ulo uhi) = ulo <  xlo && xhi <  uhi

disjoint  :: Ord r => Box r -> Box r -> Bool
{-# SPECIALIZE disjoint :: Box Double -> Box Double -> Bool #-}
{-# SPECIALIZE disjoint :: Box Rational -> Box Rational -> Bool #-}
disjoint  (Box xlo xhi) (Box ulo uhi) = xlo >= uhi || xhi <= ulo

disjoint' :: Ord r => Box r -> Box r -> Bool
{-# SPECIALIZE disjoint' :: Box Double -> Box Double -> Bool #-}
{-# SPECIALIZE disjoint' :: Box Rational -> Box Rational -> Bool #-}
disjoint' (Box xlo xhi) (Box ulo uhi) = xlo >  uhi || xhi <  ulo

symmetric :: Num r => r -> Box r
{-# SPECIALIZE symmetric :: Double -> Box Double #-}
{-# SPECIALIZE symmetric :: Rational -> Box Rational #-}
symmetric r = Box (negate r') r' where r' = abs r

singleton :: r -> Box r
{-# SPECIALIZE singleton :: Double -> Box Double #-}
{-# SPECIALIZE singleton :: Rational -> Box Rational #-}
singleton r = Box r r

union :: Ord r => Box r -> Box r -> Box r
{-# SPECIALIZE union :: Box Double -> Box Double -> Box Double #-}
{-# SPECIALIZE union :: Box Rational -> Box Rational -> Box Rational #-}
union (Box a b) (Box u v) = Box (min a u) (max b v)
