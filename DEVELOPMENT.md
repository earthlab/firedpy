# FIREDpy Development Workflow Meeting Notes (Editable)

This outlines how we will use the **biweekly FIRED development meetings** to structure progress, assign tasks, and implement repository improvements.

---

## Meeting Information
- **Frequency:** Biweekly (every 2 weeks)
   - Alternates between **working** meetings and **update** meetings
- **Duration:** 1 hour
- **Purpose:**
  - Working meetings
      - Review open pull requests
      - Review open GitHub issues
      - Assign PRs and issues
      - Update development tracker
  - Update meetings
      - Present on progress to larger team
      - Code demonstrations

---

## Project Tracking

### **FIREDpy Development Kanban Tracker**
We use **[GitHub Projects](https://github.com/earthlab/firedpy/projects)** for tracking issues and pull requests as they move through development.

### **How to Use the Kanban Board**
1. **View the board**: Go to the FIREDpy repository → **Projects** tab → **"FIRED Development Board"**
2. **Three columns represent progress stages**:
   - **To Do:** newly created or planned tasks
   - **In Progress:** actively being worked on
   - **Done:** completed tasks
3. **Link issues to the board**:
   - In each issue's sidebar → click **Projects** → select the FIRED board → choose **Add to project**
4. **Move issues between columns** as they progress

---

## Integration Branch Workflow Review

1. **Developers** work in personal feature branches (`name-dev`).
2. All changes are **PR'd into `nonprod`** (protected branch).
3. When `nonprod` is stable, **PR `nonprod` → `main`** (review by Adam/Aashish).
4. Upon merge to `main`, **tag new version** if major updates are included.

## Versioning

### **Versioning Policy**
FIREDpy follows **semantic versioning** (`MAJOR.MINOR.PATCH`):
Read more [here](https://semver.org/)

| Update Type | When to Use | Example | Notes |
|--------------|-------------|----------|-------|
| **MAJOR** | Core algorithm improvements (ex. Koppen optimization) | `2.0.0` | Update user docs and changelog |
| **MINOR** | Individual improvements (ex. cFRP integration, speed integration) | `2.1.0` | Add to changelog |
| **PATCH** | Bug fixes, refactors, or small performance updates | `2.1.1` | No API changes |

### **Release Process**
After merging major changes into `main`, create a **release tag** in GitHub:

#### **Option 1: GitHub Web Interface**
1. Go to the repository on GitHub
2. Click **"Releases"** (on the right sidebar or main navigation)
3. Click **"Create a new release"**
4. **Tag version**: Enter version number (e.g., `v2.1.0`)
5. **Release title**: Brief description (e.g., "Added cFRP integration and fire speed computation")
6. **Description**: Detailed changelog of what's new
7. Click **"Publish release"**

#### **Option 2: Command Line**
```bash
git tag -a v2.1.0 -m "Added cFRP integration and fire speed computation"
git push origin v2.1.0
```

### **Version Tracking**
- Maintain a `CHANGELOG.md` documenting all **MAJOR** and **MINOR** updates 
- Use the **[Major Change/Release Template](.github/pull_request_template_major.md)** for version-impacting changes

---

## Enhancements

### **Enhancement Requests**
Use the **[Feature Request Template](.github/ISSUE_TEMPLATE/feature_request.md)** for suggesting new functionality.

### **Current Priority Development Items**
The following are current **top development priorities** and will remain recurring agenda items until completed:

1. **Integrate Fire Speed into Algorithm**  
   - Implement fire growth rate metrics; validate against sample datasets.
   - *Link to open issues: [Search for "fire speed" issues](https://github.com/earthlab/firedpy/issues?q=is%3Aopen+is%3Aissue+fire+speed)*
   - **Lead(s):** Jerry Gamie / Danielle Losos / Nicole Hemming-Schroeder
2. **Integrate cFRP into Algorithm**  
   - Add cFRP calculation for MODIS and VIIRS AFD.
   - *Link to open issues: [Search for "cFRP" issues](https://github.com/earthlab/firedpy/issues?q=is%3Aopen+is%3Aissue+cFRP)*
   - **Lead(s):** Nate Hofford
3. **Optimize Aggregation by Köppen Region**  
   - Tune temporal/spatial aggregation using regional climate classes.
   - *Link to open issues: [Search for "Köppen" issues](https://github.com/earthlab/firedpy/issues?q=is%3Aopen+is%3Aissue+köppen)*
   - **Lead(s)**: Kyle Manley
4. **Fix Land Cover Issue**
   - Some data have NAs for land cover
   - *Link to open issues: [Search for "Land Cover" issues](https://github.com/earthlab/firedpy/issues?q=is%3Aopen+is%3Aissue+land+cover)*
   - **Lead(s)**: Aashish Mukund
---

## Structural Improvement Tasks

### **Immediate To‑Dos (Structural Improvements)**
- Organize and archive branches
- Create dev branches for team members
- Make README for code base to make development easier
- Organize existing issues

### Additional Repo Hygiene
- Add a **Versioning** section to the README (SemVer: major/minor with examples)
- Update **Development Team List** in `CONTRIBUTING.md`
- Branch protection for `nonprod` and `main` (PR-only merges)

---

## Contributing

### **Contributing Guidelines**
See **[CONTRIBUTING.md](CONTRIBUTING.md)** for guidelines on:
- Development workflows for internal team members vs external contributors
- Branch naming conventions and development workflows
- Code standards and testing requirements 
- Templates for questions, bugs, enhancements, and pull-requests

### **Development Team Roles** (need people to update with their contributions / add previous contributors)
| Name | Role | Responsibilities | Contributions |
|------|------|------------------|---------------|
| **Adam Mahood** | Project Lead | Final reviews for merges into `main` | |
| **Aashish Mukund** | Software developer | Backup reviews and algorithm development | |
| **Maxwell Cook** | Contributor | | EarthAccess fixes, ICS-209+ integration, cFRP workflow development |
| **Kyle Manly** | Contributor | Köppen Region spatial/temporal optimization | EarthAccess fixes |
| **Nicole Hemming-Schroeder** | Contributor | | Fire speed algorithm conceptualization and development |
| **Danielle Losos** | Contributor | Fire speed testing | |
| **Jerry Gammie** | Contributor | Fire speed - python conversion | Fire speed python code |
| **Nate Hofford** | Coordinator/Contributor | Meeting facilitation, cFRP integration | |

---
