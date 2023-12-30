/*
 * ------------------------------------------------------------------
 *  Author: Alaa Mroue - 2023
 *  This code is licensed under GPLv3. See LICENCE for more information.
 * ------------------------------------------------------------------
 */

#include "GlobalHandlerClass.h"

GlobalHandlerClass::GlobalHandlerClass(int globalintX, int globalintY) {
	this->globalintX = globalintX;
	this->globalintY = globalintY;

}

int GlobalHandlerClass::get_group_id(int index) {
	if (index == 0 ){
		return this->globalintX;
	}
	else if (index == 1){
		return this->globalintY;
	}

	std::cout << "get_group_id error: Index is not valid" << std::endl;
	return 0;
}
int GlobalHandlerClass::get_global_size(int index) {
	if (index == 0) {
		return this->globalintX;
	}
	else if (index == 1) {
		return this->globalintY;
	}

	std::cout << "get_global_size error: Index is not valid" << std::endl;
	return 0;
}
int GlobalHandlerClass::get_global_id(int index) {
	if (index == 0) {
		return this->globalintX;
	}
	else if (index == 1) {
		return this->globalintY;
	}

	std::cout << "get_global_id error: Index is not valid" << std::endl;
	return 0;
}
int GlobalHandlerClass::get_local_id(int index) {
	if (index == 0) {
		return this->globalintX;
	}
	else if (index == 1) {
		return this->globalintY;
	}

	std::cout << "get_local_id error: Index is not valid" << std::endl;
	return 0;
}
int GlobalHandlerClass::get_local_size(int index) {
	if (index == 0) {
		return this->globalintX;
	}
	else if (index == 1) {
		return this->globalintY;
	}

	std::cout << "get_local_size error: Index is not valid" << std::endl;
	return 0;
}