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
 * Uncopyable.h
 *
 *  Created on: 04.05.2012
 *      Author: arnem
 */

#ifndef UNCOPYABLE_H_
#define UNCOPYABLE_H_

namespace metaopt {

/**
 * creates default behavior of assignment operator that fails an assert
 */
class Uncopyable {
public:
	Uncopyable();

	virtual void operator=(Uncopyable& ref);

	virtual ~Uncopyable();
};

} /* namespace metaopt */
#endif /* UNCOPYABLE_H_ */
