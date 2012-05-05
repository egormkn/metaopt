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
