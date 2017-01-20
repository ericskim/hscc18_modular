/*
 * scots.hh
 *
 *     created: Jan 2017
 *      author: Matthias Rungger
 */

/** @file **/

#ifndef SCOTS_HH_
#define SCOTS_HH_

#include "TransitionFunction.hh"
#include "UniformGrid.hh"
#include "AbstractionGB.hh"
#include "GameSolver.hh"
#include "WinningDomain.hh"
#include "StaticController.hh"
#include "InputOutput.hh"


/* if scots is used in connection with the cudd library */
#ifdef  SCOTS_BDD
#include "dddmp.h"
#include "IndexSet.hh"
#include "SymbolicModelGB.hh"
#endif

#endif /* SCOTS_HH_ */

