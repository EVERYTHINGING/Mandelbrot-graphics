{-
    An independent implementation of the algorithm in:

    "Images of Julia sets that you can trust"
    L. H. de Figueiredo, D. Nehab, J. Stolfi, and J. B. Oliveira
    Last updated on January 8, 2013 at 10:45am.
    <http://webdoc.sub.gwdg.de/ebook/serien/e/IMPA_A/721.pdf>
-}
module Mandelbrot.Graphics.Trustworthy.Julia
  ( Color(..)
  , Image(..)
  , julia
  ) where

import Control.Monad (forM_, when)
import Control.Monad.Loops (allM, anyM)
import Data.Ix (Ix(..))
import Data.Map (Map)
import qualified Data.Map as M
import Data.Set (Set)
import qualified Data.Set as S

import Mandelbrot.Graphics.STRef.Lazy (Ref, RefM)
import qualified Mandelbrot.Graphics.STRef.Lazy as R

import Mandelbrot.Graphics.BoundingBox.D2 (Box, inside', disjoint')
import qualified Mandelbrot.Graphics.BoundingBox.D2 as B

import Mandelbrot.Graphics.QuadTree hiding (Tree, intersection)
import qualified Mandelbrot.Graphics.QuadTree as Q

data Color = Black | Gray | White | Marked
  deriving (Eq, Ord, Read, Show, Enum, Bounded, Ix)

type Cell s = Ref s Color
type Tree s = Q.Tree (Cell s)
type World r s = (Box r, Tree s, Cell s)

black :: RefM s (Tree s)
black = Leaf <$> R.new Black

gray :: RefM s (Tree s)
gray = Leaf <$> R.new Gray

white :: RefM s (Tree s)
white = Leaf <$> R.new White

-- 0
initial :: Box r -> RefM s (World r s)
{-# SPECIALIZE initial :: Box Double -> RefM s (World Double s) #-}
{-# SPECIALIZE initial :: Box Rational -> RefM s (World Rational s) #-}
initial bounds = (,,) <$> pure bounds <*> gray <*> R.new White

-- 1
refine :: Tree s -> RefM s (Tree s)
refine l@(Leaf r) = do
  c <- R.read r
  case c of
    Gray -> Quad <$> gray <*> gray <*> gray <*> gray
    _ -> pure l
refine (Quad a b c d) =
  Quad <$> refine a <*> refine b <*> refine c <*> refine d

type Graph s = (Map (Cell s) (Set (Cell s)), Map (Cell s) (Set (Cell s)))

empty :: Graph s
empty = (M.empty, M.empty)

insert :: Ref s (Graph s) -> Cell s -> Cell s -> RefM s ()
insert g a b = R.modify g $ \(as, bs) ->
  ( M.insertWith S.union a (S.singleton b) as
  , M.insertWith S.union b (S.singleton a) bs
  )

intersection :: RealFrac r => World r s -> Box r -> [Cell s]
{-# SPECIALIZE intersection :: World Double s -> Box Double -> [Cell s] #-}
{-# SPECIALIZE intersection :: World Rational s -> Box Rational -> [Cell s] #-}
intersection (bounds, tree, exterior) target
  | not (target `inside'` bounds) = [exterior] ++ Q.intersection bounds tree target
  | otherwise = Q.intersection bounds tree target

-- 2
graph :: RealFrac r => (Box r -> Box r) -> World r s -> RefM s (Graph s)
{-# SPECIALIZE graph :: (Box Double -> Box Double) -> World Double s -> RefM s (Graph s) #-}
{-# SPECIALIZE graph :: (Box Rational -> Box Rational) -> World Rational s -> RefM s (Graph s) #-}
graph f world@(bounds, tree, _) = do
  g <- R.new empty
  let go x (Leaf r) = do
        c <- R.read r
        case c of
          Gray -> forM_ (intersection world (f x)) $ insert g r
          _ -> pure ()
      go x (Quad a b c d) = do
        let (xa, xb, xc, xd) = B.split x
        go xa a
        go xb b
        go xc c
        go xd d
  go bounds tree
  R.read g

-- 3
whiten :: Graph s -> RefM s ()
whiten (fwd, rev) = mapM_ visit (M.keys fwd)
  where
    visit a = do
      c <- R.read a
      when (c == Gray) $ do
        w <- allM (\r -> (White ==) <$> R.read r) (S.toList $ M.findWithDefault S.empty a fwd)
        when w $ do
          R.write a White
          mapM_ visit (S.toList $ M.findWithDefault S.empty a rev)

-- 4
blacken :: Graph s -> RefM s ()
blacken (fwd, rev) = mapM_ visit (M.keys fwd) >> mapM_ unmark (M.keys fwd)
  where
    visit a = do
      c <- R.read a
      when (c == Gray) $ do
        w <- anyM (\r -> (`elem` [White, Marked]) <$> R.read r) (S.toList $ M.findWithDefault S.empty a fwd)
        when w $ do
          R.write a Marked
          mapM_ visit (S.toList $ M.findWithDefault S.empty a rev)
    unmark a = R.modify a u
    u White  = White
    u Gray   = Black
    u Black  = Black
    u Marked = Gray

-- 5
prune :: Tree s -> RefM s (Tree s)
prune l@(Leaf _) = pure l
prune (Quad a b c d) = do
  q <- Quad <$> prune a <*> prune b <*> prune c <*> prune d
  case q of
    Quad (Leaf a) (Leaf b) (Leaf c) (Leaf d) -> do
      a <- R.read a
      b <- R.read b
      c <- R.read c
      d <- R.read d
      case (a, b, c, d) of
        (Black, Black, Black, Black) -> black
        (White, White, White, White) -> white
        _ -> pure q
    _ -> pure q

-- 1-5
step :: RealFrac r => (Box r -> Box r) -> World r s -> RefM s (World r s)
{-# SPECIALIZE step :: (Box Double -> Box Double) -> World Double s -> RefM s (World Double s) #-}
{-# SPECIALIZE step :: (Box Rational -> Box Rational) -> World Rational s -> RefM s (World Rational s) #-}
step f (bounds, tree, exterior) = do
  tree <- refine tree
  g <- graph f (bounds, tree, exterior)
  whiten g
  blacken g
  tree <- prune tree
  pure (bounds, tree, exterior)

type Image r = Q.Tree (Box r, Color)

image :: RealFrac r => World r s -> RefM s (Image r)
{-# SPECIALIZE image :: World Double s -> RefM s (Image Double) #-}
{-# SPECIALIZE image :: World Rational s -> RefM s (Image Rational) #-}
image (bounds, tree, _) = boxed bounds <$> go tree
  where
    go (Leaf l) = Leaf <$> R.read l
    go (Quad a b c d) = Quad <$> go a <*> go b <*> go c <*> go d

-- 0-6
julia :: RealFrac r => (Box r -> Box r) -> r -> [Image r]
{-# SPECIALIZE julia :: (Box Double -> Box Double) -> Double -> [Image Double] #-}
{-# SPECIALIZE julia :: (Box Rational -> Box Rational) -> Rational -> [Image Rational] #-}
julia f r = R.run $ do
  world <- initial (B.symmetric r)
  go world
  where
    go world = do
      i <- image world
      world <- step f world
      (i :) <$> go world
