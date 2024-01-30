{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE DeriveFunctor #-}
{-# LANGUAGE DeriveFoldable #-}
{-# LANGUAGE DeriveTraversable #-}
module Mandelbrot.Graphics.QuadTree where

import Mandelbrot.Graphics.BoundingBox.D2

data Tree a = Leaf a | Quad (Tree a) (Tree a) (Tree a) (Tree a)
  deriving (Eq, Ord, Read, Show, Functor, Foldable, Traversable)

boxed :: RealFrac r => Box r -> Tree a -> Tree (Box r, a)
{-# SPECIALIZE boxed :: Box Double -> Tree a -> Tree (Box Double, a) #-}
{-# SPECIALIZE boxed :: Box Rational -> Tree a -> Tree (Box Rational, a) #-}
boxed !z (Leaf a) = Leaf (z, a)
boxed !z (Quad a b c d) =
    Quad (boxed za a) (boxed zb b) (boxed zc c) (boxed zd d)
  where
    (za, zb, zc, zd) = split z

intersection :: RealFrac r => Box r -> Tree a -> Box r -> [a]
{-# SPECIALIZE intersection :: Box Double -> Tree a -> Box Double -> [a] #-}
{-# SPECIALIZE intersection :: Box Rational -> Tree a -> Box Rational -> [a] #-}
intersection bounds tree target = go bounds tree
  where
    go z t
      | target `disjoint'` z = []
      | otherwise = case t of
          Leaf l -> [l]
          Quad a b c d ->
            let (za, zb, zc, zd) = split z
            in  concat
                  [ go za a
                  , go zb b
                  , go zc c
                  , go zd d
                  ]
