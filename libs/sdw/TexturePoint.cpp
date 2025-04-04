#include "TexturePoint.h"

TexturePoint::TexturePoint() = default;
TexturePoint::TexturePoint(float xPos, float yPos) : x(xPos), y(yPos) {}

std::ostream &operator<<(std::ostream &os, const TexturePoint &point) {
	os << "x: " << point.x << " y: " << point.y;
	return os;
}

TexturePoint TexturePoint::operator-(const TexturePoint &rh) {
	return TexturePoint{x - rh.x, y - rh.y};
}

TexturePoint TexturePoint::operator+(const TexturePoint &rh) {
	return TexturePoint{x + rh.x, y + rh.y};
}

TexturePoint TexturePoint::operator*(float rh) {
	return TexturePoint{x * rh, y * rh};
}

TexturePoint TexturePoint::operator/(float rh) {
	return TexturePoint{x/rh, y/rh};
}