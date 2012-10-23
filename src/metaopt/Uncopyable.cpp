/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    fast-tfva - efficient thermodynamic constrained flux variability analysis.
    Copyright (C) 2012  Arne Müller, arne.mueller@fu-berlin.de

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/*
 * Uncopyable.cpp
 *
 *  Created on: 04.05.2012
 *      Author: arnem
 */

#include "Uncopyable.h"
#include "assert.h"

namespace metaopt {

Uncopyable::Uncopyable() {
	// do nothing
}

void Uncopyable::operator=(Uncopyable& ref) {
	// assignment operator must not be used
	assert(false); // copying uncopyable classes may lead to memory corruption!
}

Uncopyable::~Uncopyable() {
	// do nothing
}

} /* namespace metaopt */
