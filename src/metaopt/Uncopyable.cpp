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
