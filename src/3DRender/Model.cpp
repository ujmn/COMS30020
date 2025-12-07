#include <Model.h>
#include <fstream>
#include <map>
#include <Utils.h>
#include <RenderUtils.h>
#include <RenderStrategy.h>

Model::Model(std::string objFilePath, std::string mtlFilePath)
    : objFilePath_(objFilePath)
    , mtlFilePath_(mtlFilePath)
{
    loadFile(objFilePath_, mtlFilePath_);
	auto [ver2face, face2ver] = createVertex2FaceIndex();
	verInd2faceInd_ = std::move(ver2face);
	faceIndex2verInd_ = std::move(face2ver);
	vertexNormal_ = std::move(calculateVertexNormal());
}

Model::~Model(){

}

Model::Model(const Model &rhs)
{
    objFilePath_ 		= rhs.objFilePath_;
    mtlFilePath_ 		= rhs.mtlFilePath_;
    modelTriangles_ 	= rhs.modelTriangles_;
    mirrorFaces_ 		= rhs.mirrorFaces_;
    renderStrategy_		= rhs.renderStrategy_;
    verInd2faceInd_ 	= rhs.verInd2faceInd_;
    faceIndex2verInd_ 	= rhs.faceIndex2verInd_;
    vertexNormal_ 		= rhs.vertexNormal_;
}

void Model::loadFile(std::string objFilePath, std::string mtlFilePath)
{
	std::ifstream objFile{objFilePath, std::ios::in};
	std::ifstream mtlFile{mtlFilePath, std::ios::in};
	std::vector<glm::vec3> vertics;
	if (!objFile.is_open())
    {
        std::cerr << "Failed to open file: " << objFilePath << '\n';
    }
	if (!mtlFile.is_open())
    {
        std::cerr << "Failed to open file: " << mtlFilePath << '\n';
    }

	//parse mtl
	std::map<std::string, uint32_t> str2color;
	for (std::string str; std::getline(mtlFile, str);) {
		if (!str.empty()) {
			if(str.find("newmtl") != std::string::npos) {
				auto list = split(str, ' ');	
				auto colorName = list[1];
				std::getline(mtlFile, str);
				list = split(str, ' ');	
				uint32_t color = RenderUtils::Utils::packColor(Colour{int(std::stof(list[1])*255.0), int(std::stof(list[2])*255.0), int(std::stof(list[3])*255.0)});
				str2color[colorName] = color;
			}
		}
	}

	//parse OBJ
	uint32_t color{0};
	std::string colorName{};
	for (std::string str; std::getline(objFile, str);) {
		if (!str.empty()) {
			if (str.find("usemtl") != std::string::npos) {
				std::vector<std::string> list = split(str, ' ');		
				colorName = list[1];
				color = str2color[colorName];
			}
			else if(str.front() == 'v') {
				std::vector<std::string> lineList = split(str, ' ');		
				glm::vec3 vertic{std::stof(lineList[1]), std::stof(lineList[2]), std::stof(lineList[3])};
				vertics.emplace_back(vertic);
			}else if (str.front() == 'f') {
				std::vector<std::string> lineList = split(str, ' ');		
				glm::ivec3 face{std::stoi(lineList[1]), std::stoi(lineList[2]), std::stoi(lineList[3])};
				if (colorName.empty()) {
					colorName = "Red";
				}
				color = str2color[colorName];
				auto c = RenderUtils::Utils::unpackColor(color);
				c.name = colorName;
				ModelTriangle mt{vertics[face.x-1], vertics[face.y-1], vertics[face.z-1], c};
				mt.normal = glm::normalize(glm::cross((mt.vertices[1] - mt.vertices[0]), (mt.vertices[2] - mt.vertices[0])));
				modelTriangles_.emplace_back(mt);
			}
		}
	}
}

void Model::draw(std::shared_ptr<RenderContext> context)
{
    renderStrategy_->render(*this, context);
}

const std::vector<ModelTriangle>& Model::triangles() const
{
    return modelTriangles_;
}

void Model::setMirrorFaces(std::vector<int> mirrorFaces)
{
	mirrorFaces_ = mirrorFaces;
}

std::vector<int> Model::mirrorFaceIndex() const
{
    return mirrorFaces_;
}

std::tuple<std::vector<std::vector<int>>, std::vector<std::vector<int>>> Model::createVertex2FaceIndex() {
	std::ifstream file{objFilePath_, std::ios::in};
	if (!file.is_open())
    {
        std::cerr << "Failed to open file: " << objFilePath_ << '\n';
    }

	int verIndex = 0;
	int faceIndex = 0;
	std::vector<std::vector<int>> verInd2faceInd;
	std::vector<std::vector<int>> faceIndex2verInd;
	for (std::string str; std::getline(file, str);) {
		if (!str.empty()) {
			if(str.front() == 'v') {
				verInd2faceInd.emplace_back(std::vector<int>{});
			}else if (str.front() == 'f') {
				faceIndex2verInd.emplace_back(std::vector<int>{});
				std::vector<std::string> lineList = split(str, ' ');		
				verInd2faceInd[std::stoi(lineList[1]) - 1].emplace_back(faceIndex);
				verInd2faceInd[std::stoi(lineList[2]) - 1].emplace_back(faceIndex);
				verInd2faceInd[std::stoi(lineList[3]) - 1].emplace_back(faceIndex);
				faceIndex2verInd[faceIndex].emplace_back(std::stoi(lineList[1]) - 1);
				faceIndex2verInd[faceIndex].emplace_back(std::stoi(lineList[2]) - 1);
				faceIndex2verInd[faceIndex].emplace_back(std::stoi(lineList[3]) - 1);
				faceIndex++;
			}
		}
	}

	return std::tuple<std::vector<std::vector<int>>, std::vector<std::vector<int>>>{verInd2faceInd, faceIndex2verInd};
}

std::vector<glm::vec3> Model::calculateVertexNormal()
{
	std::vector<glm::vec3> normalArr(verInd2faceInd_.size(), glm::vec3{});

	for (int i = 0; i < verInd2faceInd_.size(); i++) {
		auto faceList = verInd2faceInd_[i];
		glm::vec3 aveNormal{0.0f, 0.0f, 0.0f};
		for (const auto face : faceList) {
			aveNormal += modelTriangles_[face].normal;
		}

		aveNormal = glm::normalize(aveNormal);
		normalArr[i] = aveNormal;
	}

	return normalArr;
}

const std::vector<std::vector<int>> &Model::verInd2faceInd() const
{
	return verInd2faceInd_;
}

const std::vector<std::vector<int>> &Model::faceIndex2verInd() const
{
	return faceIndex2verInd_;
}

const std::vector<glm::vec3>& Model::vertexNormal() const
{
	return vertexNormal_;
}
