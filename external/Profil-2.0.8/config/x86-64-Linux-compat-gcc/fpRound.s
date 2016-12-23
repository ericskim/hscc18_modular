#############################################################################
#
# Basic Interval Arithmetic Subroutines Rounding Control (x86-64 Linux ELF)
# -------------------------------------------------------------------------
#
# Copyright (C) 2008 Christian Keil
#
# This file is part of PROFIL/BIAS.
#
# PROFIL/BIAS is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,
# USA.
#
# $Id $
#
#############################################################################

		.file			"fpRound.s"
		.version	"$Revision$"

		.text

		.align		4
		.globl		_BiasRoundUp
		#.type			BiasRoundUp,@function
_BiasRoundUp:
		movsd     .L1e0(%rip),%xmm0
		movsd     %xmm0,%xmm1
		addsd     .L1em30(%rip),%xmm0
		ucomisd   %xmm1,%xmm0
		ja        .LUpExit
		stmxcsr		_cw(%rip)
		movl			_cw(%rip),%eax
		andb			$0x9f,%ah
		orb				$0x40,%ah
		movl			%eax,_cw(%rip)
		ldmxcsr		_cw(%rip)
.LUpExit:
		ret
.LUp:
		#.size			_BiasRoundUp,.LUp-BiasRoundUp

		.align		4
		.globl		_BiasRoundDown
		#.type			BiasRoundDown,@function
_BiasRoundDown:
		movsd     .Lm1e0(%rip),%xmm0
		movsd     %xmm0,%xmm1
		subsd     .L1em30(%rip),%xmm0
		ucomisd   %xmm1,%xmm0
		jb        .LDownExit
		stmxcsr		_cw(%rip)
		movl			_cw(%rip),%eax
		andb			$0x9f,%ah
		orb				$0x20,%ah
		movl			%eax,_cw(%rip)
		ldmxcsr		_cw(%rip)
.LDownExit:
		ret
.LDown:
		#.size			_BiasRoundDown,.LDown-BiasRoundDown

		.align		4
		.globl		_BiasRoundNear
		#.type			BiasRoundNear,@function
_BiasRoundNear:
		movsd     .L1e0(%rip),%xmm0
		movsd     .L1e0(%rip),%xmm1
		addsd     .L1em30(%rip),%xmm0
		ucomisd   %xmm1,%xmm0
		ja        .LNearChange
		movsd     .Lm1e0(%rip),%xmm0
		movsd     .Lm1e0(%rip),%xmm1
		subsd     .L1em30(%rip),%xmm0
		ucomisd   %xmm1,%xmm0
		jb        .LNearChange
		jmp       .LNearExit
.LNearChange:
		stmxcsr		_cw(%rip)
		movl			_cw(%rip),%eax
		andb			$0x9f,%ah
		movl			%eax,_cw(%rip)
		ldmxcsr		_cw(%rip)
.LNearExit:
		ret
.LNear:
		#.size			BiasRoundNear,.LNear-BiasRoundNear

		#.local		_cw
		.comm			_cw,4,4
.L1e0:
		.double   0e1e+0
.Lm1e0:
    .double   0e-1e+0
.L1em30:
		.double   0e1e-30

		.ident		"$Id$"
