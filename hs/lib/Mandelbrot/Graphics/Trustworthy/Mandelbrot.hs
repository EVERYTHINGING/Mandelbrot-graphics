module Mandelbrot.Graphics.Trustworthy.Mandelbrot
  ( Color(..)
  , Image(..)
  , mandelbrot
  ) where


import Control.Monad (forM_)
import Data.Foldable (toList)
import Data.Maybe (mapMaybe)
import qualified Data.Set as S

import Mandelbrot.Graphics.Trustworthy.Julia
import Mandelbrot.Graphics.Trustworthy.Quadratic

import qualified Mandelbrot.Graphics.BoundingBox.D1 as D1
import Mandelbrot.Graphics.BoundingBox.D2
import Mandelbrot.Graphics.QuadTree

mandelbrot :: RealFrac r => Int -> Box r -> [Image r]
{-# SPECIALIZE mandelbrot :: Int -> Box Double -> [Image Double] #-}
{-# SPECIALIZE mandelbrot :: Int -> Box Rational -> [Image Rational] #-}
mandelbrot depth0 bounds = go depth0 (Leaf (bounds, Gray))
  where
    go depth image = image : go (depth + 1) (prune (connected depth <$> refine image))
    er = 2
    origin = symmetric 0
    core t = intersection (symmetric er) (fmap snd t) origin
    connected depth (c, Gray)
      = (,) c
      . color
      . mapMaybe (classify . core)
      . take depth
      . julia (quadratic c)
      $ er
    connected _ p = p

refine :: RealFrac r => Image r -> Image r
{-# SPECIALIZE refine :: Image Double -> Image Double #-}
{-# SPECIALIZE refine :: Image Rational -> Image Rational #-}
refine (Leaf (z, Gray)) =
  let (a, b, c, d) = split z
      l x = Leaf (x, Gray)
  in  Quad (l a) (l b) (l c) (l d)
refine l@(Leaf _) = l
refine (Quad a b c d) = Quad (refine a) (refine b) (refine c) (refine d)

prune :: Ord r => Image r -> Image r
{-# SPECIALIZE prune :: Image Double -> Image Double #-}
{-# SPECIALIZE prune :: Image Rational -> Image Rational #-}
prune l@(Leaf _) = l
prune (Quad a b c d) = case Quad (prune a) (prune b) (prune c) (prune d) of
    q@(Quad (Leaf (za, a)) (Leaf (zb, b)) (Leaf (zc, c)) (Leaf (zd, d))) ->
      case (a, b, c, d) of
        (Black, Black, Black, Black) -> Leaf (union za zd, Black)
        (White, White, White, White) -> Leaf (union za zd, White)
        _ -> q
    q -> q

classify :: [Color] -> Maybe Color
classify cs
  | S.fromList cs == S.singleton Black = Just Black
  | S.fromList cs == S.singleton White = Just White
  | otherwise = Nothing

color :: [Color] -> Color
color (White:_) = White
color (Black:_) = Black
color _ = Gray
