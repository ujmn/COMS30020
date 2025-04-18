#include "CanvasPoint.h"

CanvasPoint::CanvasPoint() :
		texturePoint(-1, -1) {}

CanvasPoint::CanvasPoint(float xPos, float yPos) :
		x(xPos),
		y(yPos),
		depth(0.0),
		brightness(1.0),
		texturePoint(-1, -1) {}

CanvasPoint::CanvasPoint(float xPos, float yPos, float pointDepth) :
		x(xPos),
		y(yPos),
		depth(pointDepth),
		brightness(1.0),
		texturePoint(-1, -1) {}

CanvasPoint::CanvasPoint(float xPos, float yPos, float pointDepth, float pointBrightness) :
		x(xPos),
		y(yPos),
		depth(pointDepth),
		brightness(pointBrightness),
		texturePoint(-1, -1) {}

std::ostream &operator<<(std::ostream &os, const CanvasPoint &point) {
	os << "(" << point.x << ", " << point.y << ", " << point.depth << ") " << point.brightness;
	return os;
}

CanvasPoint CanvasPoint::operator-(const CanvasPoint &rh) {
	return CanvasPoint{x - rh.x, y - rh.y, depth - rh.depth};
}

CanvasPoint CanvasPoint::operator*(float num) {
	return CanvasPoint{x * num, y * num, depth * num};
}

CanvasPoint CanvasPoint::operator+(const CanvasPoint &rh)
{
	return CanvasPoint{x + rh.x, y + rh.y, depth + rh.depth};
}

CanvasPoint CanvasPoint::operator/(float num) {
	return CanvasPoint{x/num, y/num, depth/num};
}