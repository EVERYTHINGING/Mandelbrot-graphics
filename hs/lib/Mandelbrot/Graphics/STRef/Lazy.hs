{-# LANGUAGE RankNTypes #-}
{-# LANGUAGE GeneralizedNewtypeDeriving #-}
module Mandelbrot.Graphics.STRef.Lazy
  ( Ref()
  , RefM()
  , run
  , new
  , read
  , write
  , modify
  , st
  )
  where

import Prelude hiding (read)

import Control.Monad.Fix (MonadFix(..))
import Control.Monad.Fail (MonadFail(..))
import Control.Monad.State (StateT, evalStateT, get, put, lift)
import Control.Monad.ST.Lazy (ST, runST)
import Data.STRef.Lazy (STRef, newSTRef, readSTRef, writeSTRef, modifySTRef)
import Numeric.Natural (Natural)

data Ref s a = Ref !Natural !(STRef s a)
instance Eq (Ref s a) where Ref _ i == Ref _ j = i == j
instance Ord (Ref s a) where Ref i _ `compare` Ref j _ = i `compare` j

newtype RefM s a = RefM{ runRefM :: StateT Natural (ST s) a }
  deriving (Functor, Applicative, Monad, MonadFix, MonadFail)

run :: (forall s . RefM s a) -> a
run a = runST (evalStateT (runRefM a) 0)

next :: RefM s Natural
next = RefM $ do
  i <- get
  put $! succ i
  pure i

new :: a -> RefM s (Ref s a)
new a = Ref <$> next <*> st (newSTRef a)

read :: Ref s a -> RefM s a
read (Ref _ r) = st (readSTRef r)

write :: Ref s a -> a -> RefM s ()
write (Ref _ r) a = st (writeSTRef r a)

modify :: Ref s a -> (a -> a) -> RefM s ()
modify (Ref _ r) f = st (modifySTRef r f)

st :: ST s a -> RefM s a
st a = RefM (lift a)
