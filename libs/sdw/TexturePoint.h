#pragma once

#include <iostream>

struct TexturePoint {
	float x{};
	float y{};

	TexturePoint();
	TexturePoint(float xPos, float yPos);
	friend std::ostream &operator<<(std::ostream &os, const TexturePoint &point);
	TexturePoint operator-(const TexturePoint &rh);
	TexturePoint operator+(const TexturePoint &rh);
	TexturePoint operator*(float rh);
	TexturePoint operator/(float rh);
};

