#include<map>
#include<string>
#include<memory>
#include<sstream>

struct SpeciesTree {
	struct Node {
		shared_ptr<Node> lc, rc;
		weak_ptr<Node> p;
		int taxon = -1;
		string name;
		map<string, double> attributes;
		map<string, double> temp;

		Node() {}

		double& operator[] (const string& name) {
			return attributes[name];
		}

		bool isLeaf() const {
			return taxon != -1;
		}

		shared_ptr<Node> parent() const {
			return p.lock();
		}

		shared_ptr<Node> leftChild() const {
			return lc;
		}

		shared_ptr<Node> rightChild() const {
			return rc;
		}

		int taxonID() const {
			return taxon;
		}

		string taxonName() const {
			return name;
		}
	};

	static string to_string(double v) {
		stringstream ss;
		ss << v;
		return ss.str();
	}

	shared_ptr<Node> r;
	map<string, double> temp;

	shared_ptr<Node> root() {
		return r;
	}

	shared_ptr<Node> addRoot() {
		r.reset(new Node());
		return r;
	}

	shared_ptr<Node> addLeft(shared_ptr<Node> node) {
		shared_ptr<Node> child(new Node());
		child->p = node;
		node->lc = child;
		return child;
	}

	shared_ptr<Node> addRight(shared_ptr<Node> node) {
		shared_ptr<Node> child(new Node());
		child->p = node;
		node->rc = child;
		return child;
	}

	string simpleTree(const string &LEN, const string& SUP) {
		string len;
		if (r->lc->attributes.count(LEN)) {
			if (r->rc->attributes.count(LEN)) len = string(":") + to_string((r->lc->attributes[LEN] + r->rc->attributes[LEN]) / 2);
			else len = string(":") + to_string(r->lc->attributes[LEN] / 2);
		}
		else {
			if (r->rc->attributes.count(LEN)) len = string(":") + to_string(r->rc->attributes[LEN] / 2);
			else len = "";
		}
		return string("(") + simpleSubtree(r->lc, LEN, SUP, true) + len + "," + simpleSubtree(r->rc, LEN, SUP, true) + len + ");";
	}

	string simpleSubtree(shared_ptr<Node> node, const string& LEN, const string& SUP, bool r = false) {
		string res;
		if (node->taxon == -1) {
			res = string("(") + simpleSubtree(node->lc, LEN, SUP) + "," + simpleSubtree(node->rc, LEN, SUP) + ")";
			if (node->attributes.count(SUP)) res += to_string(node->attributes[SUP]);
		}
		else {
			res = node->name;
		}
		if (!r && node->attributes.count(LEN)) {
			res += string(":") + to_string(node->attributes[LEN]);
		}
		return res;
	}

	string annotatedTree(const string& LEN) {
		string len;
		if (r->lc->attributes.count(LEN)) {
			if (r->rc->attributes.count(LEN)) len = string(":") + to_string((r->lc->attributes[LEN] + r->rc->attributes[LEN]) / 2);
			else len = string(":") + to_string(r->lc->attributes[LEN] / 2);
		}
		else {
			if (r->rc->attributes.count(LEN)) len = string(":") + to_string(r->rc->attributes[LEN] / 2);
			else len = "";
		}
		return string("(") + annotatedSubtree(r->lc, LEN, true) + len + "," + annotatedSubtree(r->rc, LEN, true) + len + ");";
	}

	string annotatedSubtree(shared_ptr<Node> node, const string& LEN, bool r = false) {
		string res;
		if (node->taxon == -1) {
			res = string("(") + annotatedSubtree(node->lc, LEN) + "," + annotatedSubtree(node->rc, LEN) + ")";
			res += annotations(node);
		}
		else {
			res = node->name;
		}
		if (!r && node->attributes.count(LEN)) {
			res += string(":") + to_string(node->attributes[LEN]);
		}
		return res;
	}

	string annotations(shared_ptr<Node> node) {
		string res;
		if (node->attributes.size() == 0) return res;
		bool first = true;
		res = "'[";
		for (const auto& e : node->attributes) {
			if (first) first = false;
			else res += ";";
			res += e.first + "=" + to_string(e.second);
		}
		res += "]'";
		return res;
	}
};